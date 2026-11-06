function sysStatus = evaluateLyapunovTheorems( ...
    f, x, vType, isRadUnbounded, vHasPosNearOrigin, dV, dVType, dVLabel, isDiscrete, useLaSalle)
%EVALUATELYAPUNOVTHEOREMS Validate V preconditions and apply Lyapunov theorems.
%   Checks whether the candidate Lyapunov function V satisfies the
%   hypotheses of Lyapunov's stability or instability theorems, and returns
%   the resulting system stability classification string.

    isVPosDef = ismember(vType, {'pd', 'lpd'});
    isInstabilityApplicable = vHasPosNearOrigin && ismember(dVType, {'pd', 'lpd'});

    % Verify candidate V preconditions
    if ~isVPosDef && ~isInstabilityApplicable
        if ~vHasPosNearOrigin
            error('Lyapunov:InvalidCandidateV', ...
                ['Candidate function V(x) is classified as "%s" and cannot be used: ', ...
                 'it has no strictly positive values in a neighborhood of x_bar, ', ...
                 'violating both the positive definiteness condition (lpd/pd) ', ...
                 'for stability theorems and the V(x) > 0 condition for the ', ...
                 'instability theorem.'], vType);
        else
            error('Lyapunov:InvalidCandidateV', ...
                ['Candidate function V(x) is classified as "%s" (not lpd or pd) ', ...
                 'and cannot be used for stability theorems; furthermore, %s is ', ...
                 'classified as "%s" (not lpd or pd), so the instability theorem ', ...
                 'cannot be applied either.'], vType, dVLabel, dVType);
        end
    end

    % Apply theorems
    if isVPosDef && ismember(dVType, {'nd', 'lnd'})
        if strcmp(vType, 'pd') && isRadUnbounded && strcmp(dVType, 'nd')
            sysStatus = 'globally asymptotically stable';
        else
            sysStatus = 'asymptotically stable';
        end
    elseif isVPosDef && ismember(dVType, {'nsd', 'lnsd'})
        if useLaSalle
            [localLaSalle, globalLaSalle] = checkKrasovskiiLaSalle(f, dV, x, isDiscrete);
            if localLaSalle
                if strcmp(vType, 'pd') && isRadUnbounded && strcmp(dVType, 'nsd') && globalLaSalle
                    sysStatus = 'globally asymptotically stable';
                else
                    sysStatus = 'asymptotically stable';
                end
            else
                sysStatus = 'stable';
            end
        else
            sysStatus = 'stable';
        end
    elseif isInstabilityApplicable
        sysStatus = 'unstable';
    else
        % V is valid (lpd or pd), but dV is lpsd, psd, or indefinite
        sysStatus = 'inconclusive';
    end
end

% =========================================================================
% Internal Helper for Krasovskii-LaSalle Invariance Principle
% =========================================================================

function [localLaSalle, globalLaSalle] = checkKrasovskiiLaSalle(f, dV, x, isDiscrete)
%CHECKKRASOVSKIILASALLE Test if Z = {x : dV(x) = 0} contains non-trivial trajectories.
    warnState = warning('off', 'all');
    cleanupObj = onCleanup(@() warning(warnState));

    localLaSalle = false;
    globalLaSalle = false;
    n = numel(x);
    dV = simplify(dV);

    % If dV is identically zero, every trajectory stays in Z
    if isequal(dV, sym(0))
        return;
    end

    % Extract unique non-constant irreducible factors of dV
    rawFactors = factor(dV);
    nonConstFactors = sym([]);
    for k = 1:numel(rawFactors)
        if ~isempty(symvar(rawFactors(k)))
            nonConstFactors(end + 1) = simplify(rawFactors(k)); %#ok<AGROW>
        end
    end
    uniqueFactors = unique(nonConstFactors);

    if isempty(uniqueFactors)
        uniqueFactors = dV;
    end

    localLaSalle = true;
    globalLaSalle = true;

    % A trajectory remains in Z iff it remains in the zero set of at least one factor
    for idx = 1:numel(uniqueFactors)
        p = uniqueFactors(idx);

        % Build initial zero-set equations for factor p
        gradP = jacobian(p, x);
        origin = zeros(n, 1);
        g0 = double(subs(gradP, x, origin));

        % If p(0) = 0 and gradP(0) = 0, check if p is semi-definite (sum of even powers)
        baseEqs = p;
        if abs(double(subs(p, x, origin))) < 1e-9 && norm(g0, inf) < 1e-9
            [pType, ~, ~] = classifyDefiniteness(p, x);
            if ismember(pType, {'pd', 'psd', 'lpd', 'lpsd', 'nd', 'nsd', 'lnd', 'lnsd'})
                % At any real zero of a semi-definite factor, its gradient must also vanish
                baseEqs = [p; gradP(:)];
            end
        end

        % Propagate equations along system trajectories for n steps/derivatives
        eqList = sym([]);
        for b = 1:numel(baseEqs)
            h = simplify(baseEqs(b));
            if isequal(h, sym(0))
                continue;
            end
            eqList(end + 1, 1) = (h == 0); %#ok<AGROW>
            for step = 1:n
                if ~isDiscrete
                    h = simplify(jacobian(h, x) * f);
                else
                    h = simplify(subs(h, x, f));
                end
                if ~isequal(h, sym(0))
                    eqList(end + 1, 1) = (h == 0); %#ok<AGROW>
                end
            end
        end

        if isempty(eqList)
            localLaSalle = false;
            globalLaSalle = false;
            return;
        end

        try
            sol = solve(eqList, x, 'Real', true, 'ReturnConditions', true);

            % Non-empty parameters indicate a continuous manifold of trajectories in Z
            if ~isempty(sol.parameters)
                localLaSalle = false;
                globalLaSalle = false;
                return;
            end

            firstVar = char(x(1));
            numSols = numel(sol.(firstVar));
            if numSols == 0
                continue;
            end

            solMat = zeros(numSols, n);
            for j = 1:n
                solMat(:, j) = double(sol.(char(x(j))));
            end

            rowNorms = sqrt(sum(solMat.^2, 2));
            % If any non-zero isolated invariant solution exists in R^n, global LaSalle fails
            if any(rowNorms > 1e-7)
                globalLaSalle = false;
            end
        catch
            localLaSalle = false;
            globalLaSalle = false;
            return;
        end
    end
end