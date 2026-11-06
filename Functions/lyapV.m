function lyapV(dot_V, x)
%LYAPV Classify system stability from dot_V assuming V is pd and radially unbounded.
%   LYAPV(dot_V, x) classifies the trajectory derivative dot_V (or increment
%   Delta_V) around the origin x_bar = 0 and determines system stability
%   under the assumption that the candidate Lyapunov function V(x) is
%   globally positive definite (pd) and radially unbounded.
%
%   Inputs:
%     dot_V - Symbolic scalar expression for dot_V(x) or Delta_V(x)
%     x     - Symbolic vector of state variables (n x 1)

    if nargin ~= 2
        error('lyapV:InvalidNumInputs', ...
            'Exactly 2 inputs are required: lyapV(dot_V, x).');
    end

    if ~isa(dot_V, 'sym') || ~isa(x, 'sym')
        error('lyapV:NonSymbolicInput', ...
            'Inputs dot_V and x must be symbolic (sym) objects.');
    end

    if ~isscalar(dot_V)
        error('lyapV:NonScalarInput', ...
            'Input dot_V must be a scalar symbolic expression.');
    end

    x = x(:);
    n = numel(x);

    extraVars = setdiff(symvar(dot_V), x.');
    if ~isempty(extraVars)
        error('lyapV:UnknownSymbolicVars', ...
            'Expression dot_V contains symbolic variables not listed in x.');
    end

    origin = zeros(n, 1);
    try
        dVAtOrigin = double(subs(dot_V, x, origin));
    catch
        error('lyapV:SingularAtOrigin', ...
            'Expression dot_V cannot be evaluated at the origin x_bar = 0.');
    end

    if ~isfinite(dVAtOrigin)
        error('lyapV:SingularAtOrigin', ...
            'Expression dot_V(0) evaluates to Inf or NaN.');
    end

    if abs(dVAtOrigin) > 1e-9
        error('lyapV:NonZeroAtOrigin', ...
            ['Expression dot_V(0) = %.4g, but at the equilibrium x_bar = 0 ', ...
             'dot_V(0) must be equal to 0.'], dVAtOrigin);
    end

    % Classify dot_V around the origin
    [dVType, ~, ~] = classifyDefiniteness(dot_V, x);

    % Evaluate Lyapunov theorem assuming V is 'pd' and radially unbounded
    vType = 'pd';
    isRadUnbounded = true;
    vHasPosNearOrigin = true;
    useLaSalle = false;
    isDiscreteFlag = false;
    fDummy = sym(zeros(n, 1));

    sysStatus = evaluateLyapunovTheorems( ...
        fDummy, x, vType, isRadUnbounded, vHasPosNearOrigin, ...
        dot_V, dVType, 'dot_V', isDiscreteFlag, useLaSalle);

    printLyapunovSummary('pd (assumed)', 'yes (assumed)', 'dot_V', dVType, sysStatus);
end