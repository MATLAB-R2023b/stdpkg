function indirlyap(f, x, varargin)
%INDIRLYAP Verify equilibrium stability via Lyapunov's Indirect Method (linearization).
%   INDIRLYAP linearizes a continuous-time or discrete-time nonlinear system
%   around the equilibrium point (x_bar, u_bar) and applies Lyapunov's
%   Indirect Method to classify the stability of x_bar.
%
%   Supported calling syntaxes:
%     1) Autonomous system (no external input u):
%        indirlyap(f, x, x_bar, isDiscrete)
%
%     2) System with input u (u automatically inferred from symvar(f)):
%        indirlyap(f, x, x_bar, u_bar, isDiscrete)
%
%     3) System with explicit symbolic input vector u and equilibrium input u_bar:
%        indirlyap(f, x, u, x_bar, u_bar, isDiscrete)
%        indirlyap(f, x, x_bar, u, u_bar, isDiscrete)
%
%   Inputs:
%     f          - Symbolic vector of system state equations (n x 1)
%     x          - Symbolic vector of state variables (n x 1)
%     u          - (Optional) Symbolic vector/scalar of input variables (m x 1)
%     x_bar      - Numeric vector representing the equilibrium state (n x 1)
%     u_bar      - (Optional) Numeric vector/scalar representing the equilibrium input (m x 1)
%     isDiscrete - Boolean/numeric flag (0 for continuous-time, ~= 0 for discrete-time)

    if nargin < 4 || nargin > 6
        error('indirlyap:InvalidNumInputs', ...
            ['Invalid number of inputs. Supported syntaxes:\n', ...
             '  indirlyap(f, x, x_bar, isDiscrete)\n', ...
             '  indirlyap(f, x, x_bar, u_bar, isDiscrete)\n', ...
             '  indirlyap(f, x, u, x_bar, u_bar, isDiscrete)\n', ...
             '  indirlyap(f, x, x_bar, u, u_bar, isDiscrete)']);
    end

    % 1. Validate equilibrium, shift to origin if needed, and linearize f
    [A, isDiscreteFlag] = linearizeSystemAtEquilibrium(f, x, varargin{:});

    % 2. Delegate eigenvalue classification and terminal display to indirlyapM
    indirlyapM(A, isDiscreteFlag);
end