function imaging_data=Find_XY_Position_Single_Session(imaging_data, tracking_data, tracking_fs)
% Find_XY_Position_Single_Session - find XY coordinates associated with imaging events
%
% imaging_data=Find_XY_Position_Single_Session(imaging_data, tracking_data, tracking_fs)
%   retrieves XY coordinates in input arg tracking_data which coincide with
%   an event in imaging_data. In both variables, time runs down the
%   columns. tracking_fs is the sampling frequency of the tracking data in
%   Hz.
%
%   Columns of imaging_data: 
%     Ca event time (s)
%     Cell ID
%     Ca event peak amplitude in MAD
%
%   Columns of tracking_data: 
%     time (s)
%     X-position of animal (cm)
%     Y-position of animal (cm)
%     smoothed speed of animal (cm/s)
%
%   After processing imaging_data will have three additional columns,
%   namely, the x and y positions and smoothed speed of the behavioral
%   tracking.

% expand imaging data
imaging_data(:, 4:6) = NaN;

% convert time stamps to single precision and compute for each
% imaging event to which behavioral time stamp it belongs
bin = discretize(single(imaging_data(:, 1)), single(tracking_data(:, 1)));
% if association fails for any imaging event, issue an error because this
% indicates a time mismatch issue, which could be fatal for analysis
badIx = ~isfinite(bin);
if any(badIx)
    error(sum(badIx) + " imaging events could not be time-matched with tracking events")
else    
    % assign
    imaging_data(:, 4:6)=tracking_data(bin, 2:4);
end