function printLyapunovSummary(vStr, radStr, dVLabel, dVStr, sysStr)
%PRINTLYAPUNOVSUMMARY Format and print Lyapunov analysis results to terminal.
%   PRINTLYAPUNOVSUMMARY(vStr, radStr) prints only V and radial unboundedness
%   (used by checkf).
%   PRINTLYAPUNOVSUMMARY(vStr, radStr, dVLabel, dVStr, sysStr) prints the
%   complete 4-line summary (used by dirlyap, dirlkls, and lyapV).

    fprintf('%24s: %s\n', 'V', vStr);
    fprintf('%24s: %s\n', 'V radially unbounded', radStr);

    if nargin == 5
        fprintf('%24s: %s\n', dVLabel, dVStr);
        fprintf('%24s: %s\n', 'sys', sysStr);
    end
end