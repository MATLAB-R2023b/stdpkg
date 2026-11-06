function [wType, isRadUnbounded, hasPosNearOrigin] = classifyDefiniteness(W, x)
%CLASSIFYDEFINITENESS Classify definiteness and radial unboundedness of W(x) at 0.
%   Assumes W(0) == 0. Classifies W(x) into one of:
%     'pd', 'psd', 'lpd', 'lpsd', 'nd', 'nsd', 'lnd', 'lnsd', 'indefinite'
%   and checks whether W(x) -> +infinity as ||x|| -> infinity.

    warnState = warning('off', 'all');
    cleanupObj = onCleanup(@() warning(warnState));

    x = x(:);
    n = numel(x);
    W = simplify(W);

    % Check radial unboundedness first
    isRadUnbounded = checkRadialUnboundedness(W, x);

    % Identically zero function is globally negative (and positive) semi-definite
    if isequal(W, sym(0))
        wType = 'nsd';
        hasPosNearOrigin = false;
        return;
    end

    origin = zeros(n, 1);
    gradW = jacobian(W, x);
    g0 = double(subs(gradW, x, origin));

    % Non-zero gradient at the origin implies W changes sign in every neighborhood
    if norm(g0, inf) > 1e-9
        wType = 'indefinite';
        hasPosNearOrigin = true;
        return;
    end

    % Evaluate Hessian matrix at the origin
    H0 = double(subs(hessian(W, x), x, origin));
    H0 = 0.5 * (H0 + H0.');
    [eigVecs, eigMat] = eig(H0);
    eigVals = diag(eigMat);

    eigTol = 1e-8;
    U = buildDirectionMatrix(n, eigVecs);

    % Step 1: Local classification in a ball B_r around the origin
    if all(eigVals > eigTol)
        localType = 'lpd';
        hasPosNearOrigin = true;
    elseif all(eigVals < -eigTol)
        localType = 'lnd';
        hasPosNearOrigin = false;
    elseif any(eigVals > eigTol) && any(eigVals < -eigTol)
        wType = 'indefinite';
        hasPosNearOrigin = true;
        return;
    else
        [localType, hasPosNearOrigin] = classifyDegenerateLocal(W, x, U, eigVecs, eigVals);
        if strcmp(localType, 'indefinite')
            wType = 'indefinite';
            return;
        end
    end

    % Step 2: Test if local definiteness holds globally on all R^n \ {0}
    isGlobal = verifyGlobalSign(W, x, gradW, U, localType);

    if isGlobal
        switch localType
            case 'lpd',  wType = 'pd';
            case 'lpsd', wType = 'psd';
            case 'lnd',  wType = 'nd';
            case 'lnsd', wType = 'nsd';
            otherwise,   wType = localType;
        end
    else
        wType = localType;
    end

    if ~strcmp(wType, 'pd') && ~strcmp(wType, 'psd')
        isRadUnbounded = false;
    end
end

% =========================================================================
% Internal Subfunctions for Definiteness and Radial Unboundedness
% =========================================================================

function [localType, hasPosNearOrigin] = classifyDegenerateLocal(W, x, U, eigVecs, eigVals)
%CLASSIFYDEGENERATELOCAL Analyze higher-order terms when H(0) is semidefinite or zero.
    n = numel(x);
    rSym = sym('r_pos', 'real');

    testDirs = [eigVecs, eye(n), ones(n, 1) / sqrt(n)];
    hasPos = any(eigVals > 1e-8);
    hasNeg = any(eigVals < -1e-8);
    hasZeroRay = false;

    for k = 1:size(testDirs, 2)
        uk = testDirs(:, k);
        if norm(uk) < 1e-12
            continue;
        end
        uk = uk / norm(uk);
        wRay = simplify(subs(W, x, rSym * uk));
        if isequal(wRay, sym(0))
            hasZeroRay = true;
            continue;
        end

        tSeries = taylor(wRay, rSym, 'ExpansionPoint', 0, 'Order', 14);
        if isequal(tSeries, sym(0))
            hasZeroRay = true;
            continue;
        end

        [coeffsVec, termsVec] = coeffs(tSeries, rSym);
        numTerms = numel(termsVec);
        degs = zeros(1, numTerms);
        for idx = 1:numTerms
            degs(idx) = double(polynomialDegree(termsVec(idx), rSym));
        end
        [minDeg, minIdx] = min(degs);
        leadCoeff = double(coeffsVec(minIdx));

        if mod(minDeg, 2) == 1
            localType = 'indefinite';
            hasPosNearOrigin = true;
            return;
        else
            if leadCoeff > 1e-12
                hasPos = true;
            elseif leadCoeff < -1e-12
                hasNeg = true;
            end
        end
    end

    if hasPos && hasNeg
        localType = 'indefinite';
        hasPosNearOrigin = true;
        return;
    end

    % Numerical evaluation on small concentric spheres
    smallRadii = [1e-3, 5e-3, 1e-2, 2.5e-2];
    for r = smallRadii
        vals = evalSymAtPoints(W, x, r * U);
        scale = max(abs(vals));
        tol = max(1e-13, 1e-7 * scale);
        if any(vals > tol)
            hasPos = true;
        end
        if any(vals < -tol)
            hasNeg = true;
        end
        if hasPos && hasNeg
            localType = 'indefinite';
            hasPosNearOrigin = true;
            return;
        end
    end

    % Check for curved zero-paths passing through the origin
    [hasCurvedZero, signChangesAcrossCurve] = checkLocalZeroCurves(W, x);
    if signChangesAcrossCurve
        localType = 'indefinite';
        hasPosNearOrigin = true;
        return;
    end
    if hasCurvedZero
        hasZeroRay = true;
    end

    if hasPos && ~hasNeg
        hasPosNearOrigin = true;
        if hasZeroRay
            localType = 'lpsd';
        else
            localType = 'lpd';
        end
    elseif hasNeg && ~hasPos
        hasPosNearOrigin = false;
        if hasZeroRay
            localType = 'lnsd';
        else
            localType = 'lnd';
        end
    else
        localType = 'lnsd';
        hasPosNearOrigin = false;
    end
end

function [hasZero, signChange] = checkLocalZeroCurves(W, x)
%CHECKLOCALZEROCURVES Detect non-trivial real zero curves passing near the origin.
    hasZero = false;
    signChange = false;
    n = numel(x);
    if n == 1
        return;
    end

    for i = 1:n
        try
            sol = solve(W == 0, x(i), 'Real', true);
            if isempty(sol)
                continue;
            end
            otherIdx = setdiff(1:n, i);
            otherVars = x(otherIdx);
            testPt = 1e-2 * ones(numel(otherVars), 1);

            for s = 1:numel(sol)
                val = double(subs(sol(s), otherVars, testPt));
                if isfinite(val) && isreal(val)
                    fullPt = zeros(n, 1);
                    fullPt(otherIdx) = testPt;
                    fullPt(i) = val;
                    wCheck = double(subs(W, x, fullPt));
                    if norm(fullPt) > 1e-5 && norm(fullPt) < 0.1 && abs(wCheck) < 1e-10
                        hasZero = true;
                        % Perturb around the zero curve to detect possible sign change
                        for factorPerturb = [0.5, 1.5]
                            pertPt = fullPt;
                            pertPt(i) = factorPerturb * val;
                            if norm(pertPt) > 1e-6
                                wPert = double(subs(W, x, pertPt));
                                wBase = double(subs(W, x, [testPt(1); zeros(n-1, 1)]));
                                if isfinite(wPert) && isfinite(wBase) && (wPert * wBase < -1e-16)
                                    signChange = true;
                                    return;
                                end
                            end
                        end
                    end
                end
            end
        catch
            % Continue if symbolic solve does not isolate variable i
        end
    end
end

function isGlobal = verifyGlobalSign(W, x, gradW, U, localType)
%VERIFYGLOBALSIGN Check whether local sign definiteness holds on all R^n \ {0}.
    isGlobal = true;
    n = numel(x);

    % 1. Multi-scale radial shell evaluation
    radii = unique([logspace(-2, 3, 35), 0.1, 0.25, 0.5, 0.75, 1.0, ...
                    1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 25.0, 50.0, 100.0]);

    for r = radii
        vals = evalSymAtPoints(W, x, r * U);
        scale = max(1, max(abs(vals)));
        tol = 1e-9 * scale;

        switch localType
            case 'lpd'
                if any(vals <= tol)
                    isGlobal = false;
                    return;
                end
            case 'lpsd'
                if any(vals < -tol)
                    isGlobal = false;
                    return;
                end
            case 'lnd'
                if any(vals >= -tol)
                    isGlobal = false;
                    return;
                end
            case 'lnsd'
                if any(vals > tol)
                    isGlobal = false;
                    return;
                end
        end
    end

    % 2. Symbolic check for non-zero roots along test rays
    rSym = sym('r_glob', 'real');
    for k = 1:min(size(U, 2), 4 * n + 2)
        uk = U(:, k);
        wRay = simplify(subs(W, x, rSym * uk));
        try
            rootsRay = double(solve(wRay == 0, rSym, 'Real', true));
            rootsRay = rootsRay(isfinite(rootsRay) & abs(rootsRay) > 1e-6);
            if ~isempty(rootsRay)
                if ismember(localType, {'lpd', 'lnd'})
                    isGlobal = false;
                    return;
                end
                for rt = rootsRay(:).'
                    wBeyond = double(subs(wRay, rSym, 1.1 * rt));
                    if (strcmp(localType, 'lpsd') && wBeyond < -1e-9) || ...
                       (strcmp(localType, 'lnsd') && wBeyond > 1e-9)
                        isGlobal = false;
                        return;
                    end
                end
            end
        catch
            % Ignore if ray equation is not analytically solvable
        end
    end

    % 3. Check real non-zero critical points (gradW == 0)
    try
        critSol = solve(gradW == 0, x, 'Real', true, 'ReturnConditions', true);
        if isempty(critSol.parameters)
            numPts = numel(critSol.(char(x(1))));
            for idx = 1:numPts
                pt = zeros(n, 1);
                for j = 1:n
                    pt(j) = double(critSol.(char(x(j)))(idx));
                end
                if all(isfinite(pt)) && norm(pt) > 1e-6
                    wCrit = double(subs(W, x, pt));
                    if (strcmp(localType, 'lpd') && wCrit <= 1e-9) || ...
                       (strcmp(localType, 'lpsd') && wCrit < -1e-9) || ...
                       (strcmp(localType, 'lnd') && wCrit >= -1e-9) || ...
                       (strcmp(localType, 'lnsd') && wCrit > 1e-9)
                        isGlobal = false;
                        return;
                    end
                end
            end
        end
    catch
        % Ignore if critical points cannot be solved in closed form
    end
end

function isRadUnb = checkRadialUnboundedness(V, x)
%CHECKRADIALUNBOUNDEDNESS Verify if V(x) -> +infinity as ||x|| -> infinity.
    n = numel(x);
    isRadUnb = true;

    varsInV = symvar(V);
    if numel(varsInV) < n
        isRadUnb = false;
        return;
    end

    U = buildDirectionMatrix(n, eye(n));

    % 1. Symbolic limit along primary axes and diagonals
    rSym = sym('r_lim', 'positive');
    numPrimary = min(size(U, 2), 2 * n + 4);
    for k = 1:numPrimary
        uk = U(:, k);
        vRay = subs(V, x, rSym * uk);
        try
            limVal = limit(vRay, rSym, inf);
            if ~isequal(limVal, sym(inf))
                isRadUnb = false;
                return;
            end
        catch
            % Fall back to numerical sphere growth check
        end
    end

    % 2. Numerical minimum growth check on expanding hyperspheres
    R = [10, 100, 1000, 10000];
    minVals = zeros(size(R));
    for idx = 1:numel(R)
        vals = evalSymAtPoints(V, x, R(idx) * U);
        if any(~isfinite(vals))
            isRadUnb = false;
            return;
        end
        minVals(idx) = min(vals);
    end

    if minVals(end) <= 0 || any(diff(minVals) <= 0) || (minVals(end) / minVals(end-1) < 1.05)
        isRadUnb = false;
    end
end

function U = buildDirectionMatrix(n, eigVecs)
%BUILDDIRECTIONMATRIX Construct unit direction vectors on the unit sphere in R^n.
    dirs = [eye(n), -eye(n), eigVecs, -eigVecs];

    for i = 1:n
        for j = (i + 1):n
            d1 = zeros(n, 1); d1(i) = 1; d1(j) = 1;
            d2 = zeros(n, 1); d2(i) = 1; d2(j) = -1;
            dirs = [dirs, d1, -d1, d2, -d2]; %#ok<AGROW>
        end
    end

    if n == 2
        angles = linspace(0, 2 * pi, 73);
        angles(end) = [];
        spherePts = [cos(angles); sin(angles)];
        dirs = [dirs, spherePts];
    else
        numPts = 120;
        primesList = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        spherePts = zeros(n, numPts);
        for k = 1:numPts
            pt = zeros(n, 1);
            for d = 1:n
                b = primesList(mod(d - 1, numel(primesList)) + 1);
                pt(d) = 2 * radicalInverse(k + 10, b) - 1;
            end
            if norm(pt) > 1e-6
                spherePts(:, k) = pt / norm(pt);
            else
                spherePts(1, k) = 1;
            end
        end
        dirs = [dirs, spherePts];
    end

    colNorms = sqrt(sum(dirs.^2, 1));
    validCols = colNorms > 1e-12;
    U = dirs(:, validCols) ./ colNorms(validCols);
end

function val = radicalInverse(index, base)
%RADICALINVERSE Van der Corput radical inverse for deterministic sphere sampling.
    val = 0;
    f = 1 / base;
    i = index;
    while i > 0
        val = val + f * mod(i, base);
        i = floor(i / base);
        f = f / base;
    end
end

function vals = evalSymAtPoints(expr, x, pts)
%EVALSYMATPOINTS Vectorized evaluation of symbolic scalar expr at n x M matrix pts.
    M = size(pts, 2);
    fn = matlabFunction(expr, 'Vars', num2cell(x(:).'));
    args = num2cell(pts, 2);
    vals = double(fn(args{:}));
    if isscalar(vals)
        vals = repmat(vals, 1, M);
    end
end