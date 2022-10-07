function sessionTime = getCNMFeSessionTime(filePath, ap)
% getCNMFeSessionTime(filePath, ap) - get start time of CNMFe raw files.
% 
% getCNMFeSessionTime(filePath, ap) retrieves the session start time from
% CNMFe raw data file filePath. In case the retrieval fails, NaT is
% returned.
% Note that the code relies on Inscopix software.
% Input arg ap ('analysis parameters') can be loaded via
%       [ds, ap] = getdsap;

DO_MOCKSESSIONTIMES = false;

sessionTime.start = NaT;
try
    if DO_MOCKSESSIONTIMES
        % mock code for development/testing
        warning(upper('Session times are random, generated for test purposes'))
        sessionTime.start = datetime('now');
    else
        % read session start times
        rawFile = isx.CellSet.read(char(filePath));
        sessionTime.start = rawFile.timing.start.datetime;
    end
    sessionTime.start.Format = ap.datetimeFormat;
catch ME
    formatSpec = "Determination of session times for %s failed: %s";
    warning(compose(formatSpec, filePath, ME.message))
end


