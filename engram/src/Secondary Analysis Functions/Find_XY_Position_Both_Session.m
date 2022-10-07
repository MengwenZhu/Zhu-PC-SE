function imaging_data=Find_XY_Position_Both_Session(imaging_data, tracking_data1, tracking_data2, Tcorrection, ds)
% Find_XY_Position_Both_Session - find XY coordinates associated with imaging events
%
% imaging_data=Find_XY_Position_Both_Session(imaging_data, tracking_data1, tracking_data2, tracking_fs, Tcorrection)
%   retrieves XY coordinates in input args tracking_data1 and
%   tracking_data2 which coincide with events in imaging_data.
%   tracking_data1 and tracking_data2 are assumed to be collected in
%   sessions separated by several hours of time. Time in each of these
%   variables starts at zero, and the first session is assumed to last no
%   longer than ds.Total_T_Recording seconds. imaging_data is a
%   concatenation of data from both sessions with only a single time zero,
%   namely the beginning of the first session. Events collected in the
%   second session are assumed have a time offset of at least 600 (in
%   practice, several thousand) seconds. In all variables, time runs down
%   the columns. tracking_fs is the sampling frequency of the tracking data
%   in Hz.
%
%   Columns of imaging_data: 
%     Ca event time (s)
%     Cell ID
%     Ca event peak amplitude in MAD
%
%   Columns of tracking_data1/2: 
%     time (s)
%     X-position of animal (units?)
%     Y-position of animal (units?)
%     smoothed speed of animal (units?)
%
%   After processing imaging_data will have three additional columns,
%   namely, the x and y positions and smoothed speed of the behavioral
%   tracking.


% expand imaging data
imaging_data(:, 4:6) = NaN;

% create logical index to rows corresponding to AM session
am_ix = imaging_data(:,1)<=ds.Total_T_Recording;
% then apply time correction for PM events
imaging_data(~am_ix, 1) = imaging_data(~am_ix, 1) - Tcorrection;
% use function Find_XY_Position_Single_Session to analyze both
% sessions separately
imaging_data_am = Find_XY_Position_Single_Session(imaging_data(am_ix, :), tracking_data1, ds.Noldus_Sampling_Frequency);
imaging_data_pm = Find_XY_Position_Single_Session(imaging_data(~am_ix, :), tracking_data2, ds.Noldus_Sampling_Frequency);
% unite them again, not forgetting to add Tcorrection to PM events
imaging_data_pm(:, 1) = imaging_data_pm(:, 1) + Tcorrection;
imaging_data = cat(1, imaging_data_am, imaging_data_pm);
