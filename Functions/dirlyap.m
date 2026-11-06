function dirlyap(f, x, V, x_bar, isDiscrete)
%DIRLYAP Verify system stability via Lyapunov's Direct Method (without LaSalle).
%   DIRLYAP(f, x, V, x_bar, isDiscrete) analyzes the stability of the
%   equilibrium point x_bar for the autonomous system:
%     - Continuous-time (isDiscrete == 0): x_dot = f(x)
%     - Discrete-time   (isDiscrete ~= 0): x(k+1) = f(x(k))
%   using the candidate Lyapunov function V(x), excluding the
%   Krasovskii-LaSalle invariance principle.
%
%   Inputs:
%     f          - Symbolic vector of system state equations (n x 1)
%     x          - Symbolic vector of state variables (n x 1)
%     V          - Symbolic scalar candidate Lyapunov function
%     x_bar      - Numeric vector representing the equilibrium point (n x 1)
%     isDiscrete - Boolean/numeric flag (0 for continuous-time, ~= 0 for discrete-time)

    if nargin ~= 5
        error('dirlyap:InvalidNumInputs', ...
            'Exactly 5 inputs are required: dirlyap(f, x, V, x_bar, isDiscrete).');
    end

    % 1. Validate inputs, verify equilibrium, shift to origin, and compute dV
    [fShifted, VShifted, dV, dVLabel, isDiscreteFlag] = ...
        validateAndShiftSystem(f, x, V, x_bar, isDiscrete);

    % 2. Classify candidate function V and trajectory change dV
    [vType, isRadUnbounded, vHasPosNearOrigin] = classifyDefiniteness(VShifted, x);
    [dVType, ~, ~] = classifyDefiniteness(dV, x);

    % 3. Apply Lyapunov Direct Method theorems (useLaSalle = false)
    useLaSalle = false;
    sysStatus = evaluateLyapunovTheorems( ...
        fShifted, x, vType, isRadUnbounded, vHasPosNearOrigin, ...
        dV, dVType, dVLabel, isDiscreteFlag, useLaSalle);

    % 4. Display results in terminal
    if isRadUnbounded
        radStr = 'yes';
    else
        radStr = 'no';
    end
    printLyapunovSummary(vType, radStr, dVLabel, dVType, sysStatus);
end