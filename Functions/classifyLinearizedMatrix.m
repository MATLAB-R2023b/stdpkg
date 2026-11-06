function sysStatus = classifyLinearizedMatrix(A, isDiscrete)
%CLASSIFYLINEARIZEDMATRIX Apply Lyapunov's Indirect Method to state matrix A.
%   Checks that A is a valid finite real square matrix and classifies the
%   equilibrium point x_bar based on the eigenvalues of A for continuous-time
%   (isDiscrete == 0) or discrete-time (isDiscrete ~= 0) systems.

    if isempty(A) || ~ismatrix(A) || size(A, 1) ~= size(A, 2)
        error('Lyapunov:InvalidMatrixA', ...
            'Input matrix A must be a non-empty square 2D matrix (n x n).');
    end

    % Allow constant symbolic matrices by converting them to double
    if isa(A, 'sym')
        if ~isempty(symvar(A))
            error('Lyapunov:SymbolicMatrixA', ...
                ['Matrix A contains free symbolic variables. Evaluate A at ', ...
                 'the equilibrium point before calling indirlyapM.']);
        end
        try
            A = double(A);
        catch
            error('Lyapunov:InvalidMatrixA', ...
                'Symbolic matrix A could not be converted to a finite numeric matrix.');
        end
    end

    if ~isnumeric(A) || ~isreal(A) || any(~isfinite(A(:)))
        error('Lyapunov:InvalidMatrixA', ...
            'Matrix A must be a finite real numeric square matrix.');
    end

    if ~isscalar(isDiscrete) || (~isnumeric(isDiscrete) && ~islogical(isDiscrete)) || ...
       ~isfinite(double(isDiscrete))
        error('Lyapunov:InvalidTimeFlag', ...
            'The final argument isDiscrete must be a finite numeric or logical scalar.');
    end
    isDiscreteFlag = (double(isDiscrete) ~= 0);

    % Compute eigenvalues of A
    eigVals = eig(double(A));
    eigTol = 1e-9;

    if ~isDiscreteFlag
        % Continuous-time Indirect Method:
        % - Re(lambda_i) < 0 for all i         -> asymptotically stable (A is Hurwitz)
        % - Re(lambda_i) > 0 for at least one  -> unstable
        % - Re(lambda_i) <= 0 with max Re == 0 -> inconclusive
        realParts = real(eigVals);
        if all(realParts < -eigTol)
            sysStatus = 'asymptotically stable';
        elseif any(realParts > eigTol)
            sysStatus = 'unstable';
        else
            sysStatus = 'inconclusive';
        end
    else
        % Discrete-time Indirect Method:
        % - |lambda_i| < 1 for all i         -> asymptotically stable (A is Schur)
        % - |lambda_i| > 1 for at least one  -> unstable
        % - |lambda_i| <= 1 with max |.| == 1 -> inconclusive
        moduli = abs(eigVals);
        if all(moduli < 1 - eigTol)
            sysStatus = 'asymptotically stable';
        elseif any(moduli > 1 + eigTol)
            sysStatus = 'unstable';
        else
            sysStatus = 'inconclusive';
        end
    end
end