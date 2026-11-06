function [fShifted, VShifted, dV, dVLabel, isDiscreteFlag] = ...
    validateAndShiftSystem(f, x, V, x_bar, isDiscrete)
%VALIDATEANDSHIFTSYSTEM Validate inputs, check equilibrium, shift to 0, and compute dV.
%   Checks dimensions and symbolic consistency of (f, x, V, x_bar, isDiscrete),
%   verifies that x_bar is an equilibrium point of f, builds the shifted
%   auxiliary system centered at the origin if x_bar ~= 0, verifies that
%   V is class C^(1) with V(x_bar) = 0, and computes dot_V or Delta_V.

    if ~isa(f, 'sym') || ~isa(x, 'sym') || ~isa(V, 'sym')
        error('Lyapunov:NonSymbolicInput', ...
            'Inputs f, x, and V must be symbolic (sym) objects.');
    end

    x = x(:);
    f = f(:);
    x_bar = double(x_bar(:));
    n = numel(x);

    if numel(f) ~= n
        error('Lyapunov:DimensionMismatch', ...
            'The dimension of system vector f (%d) must match state vector x (%d).', ...
            numel(f), n);
    end

    if ~isscalar(V)
        error('Lyapunov:NonScalarV', ...
            'The candidate Lyapunov function V must be a scalar symbolic expression.');
    end

    if numel(x_bar) ~= n || any(~isfinite(x_bar))
        error('Lyapunov:InvalidEquilibriumVector', ...
            'x_bar must be a finite numeric vector of length %d.', n);
    end

    if ~isscalar(isDiscrete) || (~isnumeric(isDiscrete) && ~islogical(isDiscrete))
        error('Lyapunov:InvalidTimeFlag', ...
            'The final argument isDiscrete must be a numeric or logical scalar.');
    end
    isDiscreteFlag = (double(isDiscrete) ~= 0);

    extraVars = setdiff([symvar(f), symvar(V)], x.');
    if ~isempty(extraVars)
        error('Lyapunov:UnknownSymbolicVars', ...
            'Expressions f or V contain symbolic variables not listed in x.');
    end

    % Verify that x_bar is a valid equilibrium point
    try
        fAtEq = double(subs(f, x, x_bar));
    catch
        error('Lyapunov:SingularAtEquilibrium', ...
            'System function f(x) cannot be evaluated at x_bar (singular or undefined).');
    end

    if any(~isfinite(fAtEq))
        error('Lyapunov:SingularAtEquilibrium', ...
            'System function f(x_bar) evaluates to Inf or NaN.');
    end

    eqTol = 1e-9;
    if ~isDiscreteFlag
        if norm(fAtEq, inf) > eqTol
            error('Lyapunov:NotAnEquilibrium', ...
                ['The provided x_bar is not an equilibrium point for the ', ...
                 'continuous-time system: f(x_bar) ~= 0 (max residual = %.3e).'], ...
                norm(fAtEq, inf));
        end
    else
        if norm(fAtEq - x_bar, inf) > eqTol
            error('Lyapunov:NotAnEquilibrium', ...
                ['The provided x_bar is not an equilibrium point for the ', ...
                 'discrete-time system: f(x_bar) ~= x_bar (max residual = %.3e).'], ...
                norm(fAtEq - x_bar, inf));
        end
    end

    % Build auxiliary system translating the equilibrium to the origin if x_bar ~= 0
    if any(abs(x_bar) > eqTol)
        if ~isDiscreteFlag
            fShifted = simplify(subs(f, x, x + x_bar));
        else
            fShifted = simplify(subs(f, x, x + x_bar) - x_bar);
        end
        VShifted = simplify(subs(V, x, x + x_bar));
    else
        fShifted = simplify(f);
        VShifted = simplify(V);
    end

    origin = zeros(n, 1);

    % Check C^(1) regularity of V around the origin
    vStr = char(VShifted);
    nonSmoothTokens = {'abs(', 'sign(', 'heaviside(', 'dirac(', 'floor(', 'ceil(', 'round('};
    for k = 1:numel(nonSmoothTokens)
        if contains(vStr, nonSmoothTokens{k})
            error('Lyapunov:NotClassC1', ...
                ['Candidate function V(x) is not usable: it contains ', ...
                 'non-differentiable term "%s" violating C^(1) regularity.'], ...
                nonSmoothTokens{k});
        end
    end

    try
        vAtOrigin = double(subs(VShifted, x, origin));
        gradV = jacobian(VShifted, x);
        gAtOrigin = double(subs(gradV, x, origin));
    catch
        error('Lyapunov:NotClassC1', ...
            ['Candidate function V(x) is not usable: V(x) or its gradient ', ...
             'cannot be evaluated at x_bar (violating C^(1) regularity).']);
    end

    if ~isfinite(vAtOrigin) || any(~isfinite(gAtOrigin))
        error('Lyapunov:NotClassC1', ...
            ['Candidate function V(x) is not usable: V(x) or its gradient ', ...
             'evaluates to Inf or NaN at x_bar (violating C^(1) regularity).']);
    end

    if abs(vAtOrigin) > eqTol
        error('Lyapunov:NonZeroAtEquilibrium', ...
            ['Candidate function V(x) is not usable: V(x_bar) = %.4g, ', ...
             'but Lyapunov theorems require V(x_bar) = 0.'], vAtOrigin);
    end

    % Compute trajectory derivative (continuous-time) or increment (discrete-time)
    if ~isDiscreteFlag
        dV = simplify(gradV * fShifted);
        dVLabel = 'dot_V';
    else
        dV = simplify(subs(VShifted, x, fShifted) - VShifted);
        dVLabel = 'Delta_V';
    end
end