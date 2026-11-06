function indirlyapM(A, isDiscrete)
%INDIRLYAPM Verify equilibrium stability via Lyapunov's Indirect Method from matrix A.
%   INDIRLYAPM(A, isDiscrete) determines the stability of the equilibrium
%   point x_bar of a nonlinear system given its linearized state matrix A:
%     - Continuous-time (isDiscrete == 0): Hurwitz criterion on eig(A)
%     - Discrete-time   (isDiscrete ~= 0): Schur criterion on eig(A)
%
%   Inputs:
%     A          - Numeric (or constant symbolic) square matrix (n x n)
%     isDiscrete - Boolean/numeric flag (0 for continuous-time, ~= 0 for discrete-time)

    if nargin ~= 2
        error('indirlyapM:InvalidNumInputs', ...
            'Exactly 2 inputs are required: indirlyapM(A, isDiscrete).');
    end

    % Classify stability of x_bar from the eigenvalues of A
    sysStatus = classifyLinearizedMatrix(A, isDiscrete);

    % Print formatted terminal output
    fprintf('%24s: %s\n', 'x_bar', sysStatus);
end