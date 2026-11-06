function [A, isDiscreteFlag] = linearizeSystemAtEquilibrium(f, x, varargin)
%LINEARIZESYSTEMATEQUILIBRIUM Validate equilibrium (x_bar, u_bar) and linearize f.
%   Parses optional input arguments (u, u_bar), verifies that (x_bar, u_bar)
%   satisfies the continuous-time (f(x_bar, u_bar) == 0) or discrete-time
%   (f(x_bar, u_bar) == x_bar) equilibrium condition, constructs the shifted
%   auxiliary system centered at the origin if x_bar ~= 0 or u_bar ~= 0, and
%   computes the linearized state matrix A = df/dx evaluated at the equilibrium.

    if ~isa(f, 'sym') || ~isa(x, 'sym')
        error('Lyapunov:NonSymbolicInput', ...
            'System vector f and state vector x must be symbolic (sym) objects.');
    end

    x = x(:);
    f = f(:);
    n = numel(x);

    if numel(f) ~= n
        error('Lyapunov:DimensionMismatch', ...
            'The dimension of system vector f (%d) must match state vector x (%d).', ...
            numel(f), n);
    end

    % Parse variable argument list: (x_bar, isDiscrete), (x_bar, u_bar, isDiscrete),
    % (u, x_bar, u_bar, isDiscrete), or (x_bar, u, u_bar, isDiscrete)
    numExtra = numel(varargin);
    u = sym([]);
    u_bar = [];

    if numExtra == 2
        % Syntax: (f, x, x_bar, isDiscrete)
        x_bar = varargin{1};
        isDiscrete = varargin{2};
    elseif numExtra == 3
        % Syntax: (f, x, x_bar, u_bar, isDiscrete)
        x_bar = varargin{1};
        u_bar = varargin{2};
        isDiscrete = varargin{3};
        u = setdiff(symvar(f), x.');
        u = u(:);
    elseif numExtra == 4
        if isa(varargin{1}, 'sym')
            % Syntax: (f, x, u, x_bar, u_bar, isDiscrete)
            u = varargin{1}(:);
            x_bar = varargin{2};
            u_bar = varargin{3};
            isDiscrete = varargin{4};
        elseif isa(varargin{2}, 'sym')
            % Syntax: (f, x, x_bar, u, u_bar, isDiscrete)
            x_bar = varargin{1};
            u = varargin{2}(:);
            u_bar = varargin{3};
            isDiscrete = varargin{4};
        else
            error('Lyapunov:InvalidInputSyntax', ...
                'When calling indirlyap with 6 inputs, u must be a symbolic variable/vector.');
        end
    end

    % Validate x_bar
    if ~isnumeric(x_bar) || isempty(x_bar)
        error('Lyapunov:InvalidEquilibriumVector', ...
            'Equilibrium state x_bar must be a numeric vector of length %d.', n);
    end
    x_bar = double(x_bar(:));
    if numel(x_bar) ~= n || ~isreal(x_bar) || any(~isfinite(x_bar))
        error('Lyapunov:InvalidEquilibriumVector', ...
            'x_bar must be a finite real numeric vector of length %d.', n);
    end

    % Validate isDiscrete flag
    if ~isscalar(isDiscrete) || (~isnumeric(isDiscrete) && ~islogical(isDiscrete)) || ...
       ~isfinite(double(isDiscrete))
        error('Lyapunov:InvalidTimeFlag', ...
            'The final argument isDiscrete must be a finite numeric or logical scalar.');
    end
    isDiscreteFlag = (double(isDiscrete) ~= 0);

    % Validate u and u_bar consistency
    remainingVars = setdiff(symvar(f), [x(:).', u(:).']);
    if ~isempty(remainingVars)
        error('Lyapunov:MissingInputEquilibrium', ...
            ['System f contains symbolic variables (%s) not included in x. ', ...
             'Please pass u_bar (and optionally u) to specify the equilibrium input.'], ...
            char(remainingVars(1)));
    end

    if ~isempty(u_bar)
        if ~isnumeric(u_bar) || ~isreal(u_bar) || any(~isfinite(u_bar(:)))
            error('Lyapunov:InvalidEquilibriumInput', ...
                'Equilibrium input u_bar must be a finite real numeric vector or scalar.');
        end
        u_bar = double(u_bar(:));

        if isempty(u)
            % Function f does not explicitly depend on u, accept scalar/vector u_bar
            u = sym('u_dummy', [numel(u_bar), 1], 'real');
        elseif numel(u) ~= numel(u_bar)
            error('Lyapunov:InputDimensionMismatch', ...
                'The length of u_bar (%d) must match the number of input variables u (%d).', ...
                numel(u_bar), numel(u));
        end
    elseif ~isempty(u)
        error('Lyapunov:MissingInputEquilibrium', ...
            'Symbolic input u was provided, but numeric equilibrium input u_bar is missing.');
    end

    % Check for non-differentiable functions violating C^(1) hypothesis
    fStr = char(f);
    nonSmoothTokens = {'abs(', 'sign(', 'heaviside(', 'dirac(', 'floor(', 'ceil(', 'round('};
    for k = 1:numel(nonSmoothTokens)
        if contains(fStr, nonSmoothTokens{k})
            error('Lyapunov:NotClassC1', ...
                ['System f contains non-differentiable term "%s" violating ', ...
                 'the C^(1) differentiability requirement for linearization.'], ...
                nonSmoothTokens{k});
        end
    end

    % Substitute equilibrium input u = u_bar and check equilibrium condition at x = x_bar
    allVars = [x; u];
    allEqVals = [x_bar; u_bar];

    try
        fAtEq = double(subs(f, allVars, allEqVals));
    catch
        error('Lyapunov:SingularAtEquilibrium', ...
            'System function f cannot be evaluated at (x_bar, u_bar) (division by zero or singularity).');
    end

    if any(~isfinite(fAtEq)) || ~isreal(fAtEq)
        error('Lyapunov:SingularAtEquilibrium', ...
            'System function f(x_bar, u_bar) evaluates to Inf, NaN, or complex values.');
    end

    eqTol = 1e-9;
    if ~isDiscreteFlag
        % Continuous-time equilibrium condition: f(x_bar, u_bar) = 0
        if norm(fAtEq, inf) > eqTol
            error('Lyapunov:NotAnEquilibrium', ...
                ['The pair (x_bar, u_bar) is not an equilibrium for the ', ...
                 'continuous-time system: f(x_bar, u_bar) ~= 0 (max residual = %.3e).'], ...
                norm(fAtEq, inf));
        end
    else
        % Discrete-time equilibrium condition: f(x_bar, u_bar) = x_bar
        if norm(fAtEq - x_bar, inf) > eqTol
            error('Lyapunov:NotAnEquilibrium', ...
                ['The pair (x_bar, u_bar) is not an equilibrium for the ', ...
                 'discrete-time system: f(x_bar, u_bar) ~= x_bar (max residual = %.3e).'], ...
                norm(fAtEq - x_bar, inf));
        end
    end

    % Construct auxiliary system shifted to the origin if x_bar ~= 0 or u_bar ~= 0
    if ~isempty(u)
        fSubU = subs(f, u, u_bar);
    else
        fSubU = f;
    end

    if any(abs(x_bar) > eqTol)
        if ~isDiscreteFlag
            fShifted = subs(fSubU, x, x + x_bar);
        else
            fShifted = subs(fSubU, x, x + x_bar) - x_bar;
        end
    else
        fShifted = fSubU;
    end

    % Compute Jacobian matrix A around the origin of the shifted system
    origin = zeros(n, 1);
    try
        symJac = jacobian(fShifted, x);
        A = double(subs(symJac, x, origin));
    catch
        error('Lyapunov:NotDifferentiableAtEquilibrium', ...
            'The Jacobian df/dx cannot be evaluated at (x_bar, u_bar) (not C^(1) at equilibrium).');
    end

    if any(~isfinite(A(:))) || ~isreal(A)
        error('Lyapunov:NotDifferentiableAtEquilibrium', ...
            'The Jacobian matrix A at (x_bar, u_bar) contains Inf, NaN, or complex entries.');
    end
end