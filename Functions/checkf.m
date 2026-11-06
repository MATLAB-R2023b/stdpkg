function checkf(V, x, x_bar)
%CHECKF Classify definiteness and radial unboundedness of a scalar function.
%   CHECKF(V, x, x_bar) evaluates whether the scalar symbolic function V(x)
%   is (locally or globally) positive/negative (semi-)definite or indefinite
%   around the evaluation point x_bar, and checks if it is radially unbounded.
%
%   Inputs:
%     V     - Symbolic scalar function to classify (e.g., V, dot_V, Delta_V)
%     x     - Symbolic vector of state variables (n x 1)
%     x_bar - Numeric vector representing the evaluation point (n x 1)

    if nargin ~= 3
        error('checkf:InvalidNumInputs', ...
            'Exactly 3 inputs are required: checkf(V, x, x_bar).');
    end

    if ~isa(V, 'sym') || ~isa(x, 'sym')
        error('checkf:NonSymbolicInput', ...
            'Inputs V and x must be symbolic (sym) objects.');
    end

    if ~isscalar(V)
        error('checkf:NonScalarFunction', ...
            'Input function V must be a scalar symbolic expression.');
    end

    x = x(:);
    x_bar = double(x_bar(:));
    n = numel(x);

    if numel(x_bar) ~= n || any(~isfinite(x_bar))
        error('checkf:InvalidPointVector', ...
            'x_bar must be a finite numeric vector of length %d matching x.', n);
    end

    extraVars = setdiff(symvar(V), x.');
    if ~isempty(extraVars)
        error('checkf:UnknownSymbolicVars', ...
            'Function V contains symbolic variables not listed in x.');
    end

    % Shift coordinate frame if x_bar is not the origin
    eqTol = 1e-9;
    if any(abs(x_bar) > eqTol)
        VShifted = simplify(subs(V, x, x + x_bar));
    else
        VShifted = simplify(V);
    end

    origin = zeros(n, 1);
    try
        vAtOrigin = double(subs(VShifted, x, origin));
    catch
        error('checkf:NotContinuousAtPoint', ...
            'Function V(x) cannot be evaluated at x_bar (singular or undefined).');
    end

    if ~isfinite(vAtOrigin)
        error('checkf:NotContinuousAtPoint', ...
            'Function V(x_bar) evaluates to Inf or NaN.');
    end

    if abs(vAtOrigin) > eqTol
        error('checkf:NonZeroAtPoint', ...
            ['Function V(x) evaluates to %.4g at x_bar, but definiteness ', ...
             'definitions require V(x_bar) = 0.'], vAtOrigin);
    end

    % Classify definiteness and radial unboundedness
    [vType, isRadUnbounded, ~] = classifyDefiniteness(VShifted, x);

    if isRadUnbounded
        radStr = 'yes';
    else
        radStr = 'no';
    end

    printLyapunovSummary(vType, radStr);
end