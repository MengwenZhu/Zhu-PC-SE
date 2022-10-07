function [ds, ap] = getdsap
% GETDSAP - define parameters describing analysis parameters and data properties.
% [ds, ap] = getdsap outputs structures ap and ds:
%   - ds ('data set') describes the data to be used for engram analysis and 
%   associated tasks
%   - ap ('analysis parameters') govern the analysis of the data
% 
% The function needs to be called by the top-level scripts running the
% engram analysis and associated tasks.

%% Section 1: define metadata excel sheets
% 1. ds - parameters describing the data set 
% This primary analysis metadata file does not contain extracted start times of recordings
% This is only used for primary batch analysis because we haven't got to
% the step where we generate denoised_cellset files within CNMFe folders,
% which are the origin of start times extracted. Also, primary batch
% analysis does not resuire any start times to proceed.
ds.metadataFileName_forPrimaryAnalysis = 'a5-i-KO Mice Ca-Imaging MetaData For Primary Analysis.xlsx';

% This secondary analysis metadata file contains extracted start times of
% recordings (four more columns appended), which will be used later on to correct time values of
% exported calcium events in secondary analysis
ds.metadataFileName_forSecondaryAnalysis = 'a5-i-KO Mice Ca-Imaging MetaData For Secondary Analysis.xlsx';

%% Section 2: define directories to data files & machine info
% to define directories, find out on which machine we are 
% Note: please add your own machine into the list in order for the program
% to run correctly. You can just copy and paste from case 'machine name'
% until the start of another case 'machine name', and make changes
% accordingly. This is for the convenience if you work on multiple machines
% and want to quickly get started. You could know your machine name after
% you run startup.m
if isunix
    [~, machine] = system('hostname');
    machine=deblank(machine);
else
    machine=getenv('computername');
end
switch lower(machine)
    case 'anespl19'
        % directory harboring excel file containing metadata about animals
        ds.metadataPath='C:\Users\pearcelabuser\Desktop\Raw Data';
        % base directory harboring raw Ca-imaging & behavioral-tracking data
        ds.CNMFe_Inscopix_rawdataPath = 'C:\Users\pearcelabuser\Desktop\Raw Data';
    case 'anespl14'
        % directory harboring excel file containing metadata about animals
        ds.metadataPath='D:\Meta_AnaOutput_Graphs\GAD Mice Ca-Imaging MetaData';
        % base directory harboring raw Ca-imaging & behavioral-tracking data
        ds.CNMFe_Inscopix_rawdataPath = 'H:/CNMFe Analysis/CNMFe RawData Files';
    otherwise
        error('machine not defined')
end

%% Section 3: define analysis parameters

% duration of one experiment in sec
ds.Total_T_Recording = 600;

% XYposition sampling rate of Noldus (behavioral tracking) in Hz
ds.Noldus_Sampling_Frequency = 25;

% 2. ap - analysis parameters 
% execution mode
ap.exec_type = "parallel";

% time format
ap.datetimeFormat = 'yyyy-MM-dd, HH:mm:ss.SSS';

% define under which amplitude will the event be excluded from place cell analysis
ap.Amp_Exclusion_Threshold = 1;

% define the amplitude correction factor for place cell analysis
ap.Amp_factor_placeness = 1.5;

% define the amplitude correction factor for PV analysis
ap.Amp_factor_PV = 0.5;

% define what portion of Inscopix & Noldus data to include in engram
% analysis (measured in sec, e.g. last 300sec)
ap.Engram_Start_T = 0;
ap.Engram_End_T = 600;

% parameters used to define upper & lower bounds of calcium event rate
% distribution, as well as how you bin the distribution (in spikes/sec)
ap.Ca_Rate_Upper_Bound=0.5; 
ap.Ca_Rate_Lower_Bound=0; 
ap.Ca_Rate_Interval=0.02; 

% triangular filter size used to smooth Noldus behavioral tracking data in
% sec
ap.Behavior_Filter_Size = 0.5;

% the lower speed threshold in cm/sec that defines whether mouse is mobile or
% immobile, and the upper spped threshold in cm/sec that removes some
% potential behavioral tracking artifacts
ap.Low_Mobility_Threshold=2;
ap.High_Mobility_Threshold=100;

% the threshold of mobility below which a recording is not considered worth including
% in statistical/graphical results (in % of time animal spent moving)
ap.Mobility_pass_threshold = 0.15;

% the additional edges added to the matrix used to fully depict the experimental 
% arena (make sure that all detected mouse activities are included)
ap.Arena_Add_Edge=2;

% specify how to bin the arena in each dimension for PV, calcium rate map
% correlation analysis
ap.Map_Division_For_PV=14; 

% specify how to bin the arena in each dimension for place cell analysis (MI shuffle test)
ap.Map_Division_For_MI=9; 

% For MI shuffle test ONLY: number of times that Ca event timestamps will be shuffled 
% to generate a mutual information 'null' distribution)
ap.Num_MI_Shuffle = 1;

% parameters used to define upper & lower bounds of mutual information
% values distribution, as well as how you bin the distribution
ap.MI_Upper_Bound=0.02; %manually define upper bound of MI distribution
ap.MI_Lower_Bound=0; %manually define lower bound of MI distribution
ap.MI_Interval=0.001; %manually define how to bin the MI distribution

% mutual information shuffle test p-value threshold
ap.MI_p_value_Threshold=0.05;

% parameters used to define upper & lower bounds of mutual information
% shuffle test p-values distribution, as well as how you bin the distribution
ap.MI_pvalue_Upper_Bound=1; %manually define upper bound of MI p-value distribution
ap.MI_pvalue_Lower_Bound=0; %manually define lower bound of MI p-value distribution
ap.MI_pvalue_Interval=0.05; %manually define how to bin the MI p-value distribution

% define Gaussian filter used to smooth Ca rate maps that will used for
% correlation analysis as quantitative evaluation of spatial memory
ap.Gaussian_Filter=Gaussian_filter(5, 3); % first number is filter size and second is std (sigma)

% number of times you wish to shuffle calcium rate map correlation
% distribution
ap.Ca_Rate_Map_Num_Shuffle=1;

% the lowest number of place cells (defined by MI shuffle test) that are 
% active during both sessions that we choose perform correlation analysis;
% if the actual place cells are less than this threshold, we will include
% all cells for that pair of recordings.
ap.Num_Place_Cell_Threshold=10000;

% the lowest value of calcium rate map correlation that we still consider
% as memory
ap.Ca_Rate_Map_Corr_Threshold=0;

% parameters used to define upper & lower bounds of calcium rate map
% correlation distribution, as well as how you bin the distribution
ap.Ca_Rate_Map_Corr_Upper_Bound=1; 
ap.Ca_Rate_Map_Corr_Lower_Bound=-0.3; 
ap.Ca_Rate_Map_Corr_Interval=0.05;

% number of times you wish to shuffle PV correlation distribution
ap.PV_Corr_Num_Shuffle=1; % define shuffling times

% parameters used to define upper & lower bounds of PV correlation distribution, 
% as well as how you bin the distribution
ap.PV_Corr_Upper_Bound=1; 
ap.PV_Corr_Lower_Bound=-1; 
ap.PV_Corr_Interval=0.05;

% define how to bin smoothed calcium rate map that is used for graphing
% purposes (divide into how many bins in each dimension)
ap.Ca_Map_Graph_Division=49;

% define Gaussian filter used to smooth calcium rate maps that are used for
% GRAPHING ONLY (does not affect correlation memory analysis)
ap.Gaussian_Filter_For_Graph=Gaussian_filter(16, 4); % first number is filter size and second is std (sigma)

% p-value threshold from mutual information test that select which cells
% you wish to visualize when plotting place field
ap.Place_Cell_Visualization_Threshold=1;

% the lowest value of PV correlation that we still consider as memory trace
ap.PV_Corr_Threshold=0;

end