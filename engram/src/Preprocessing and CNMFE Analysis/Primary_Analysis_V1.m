function Primary_Analysis_V1(Info_Table)
%% This function uses MATLAB API in IDPS v1.6.0 to perform automated primary analysis 
% the final objective is to generate denoised-cnmfe-cellset that is used
% for longitudinal registration (LR) and event detection, which are to
% prepare for secondary analysis

% the output of the function is a bunch of files (not any variables)

% input is a row of meta-data excel sheet (Info_Table)

% These isx.function could be found in the following directory on ANESPL14: G:\Inscopix\Data Processing\+isx

% If you wonder anything about the details of these isx.function, just type
% in help isx.function within the MATLAB command window(e.g. help isx.cnmfe)
%% Step 1: Preparation
%% 1.1 set current working folder

% define the folder where calcium imaging recording could be found
project_folder = convertStringsToChars(Info_Table.Recording_Directory);
% set the folder as current folder
cd(project_folder)

% generate CNMFe directory if it does not exist yet
if ~exist(string(Info_Table.CNMFe_Directory), 'dir')
    mkdir(string(Info_Table.CNMFe_Directory));
end

%% 1.2 load the calcium imaging recording

raw_movie_file = Info_Table.Recording_Name;

%% Step 2: Preprocess the recording 
%% 2.1: make filenames that will be used later on in the analysis

pp_file = 'prepro.isxd';
bp_file = 'bp_filter.isxd';
mc_ref = 'mc_ref.isxd';
mc_file = 'motion-corr.isxd';
dff_file = 'dff.isxd';
maxdff_file = 'max_dff.isxd';
cnmfe_file = 'cnmfe-cellset.isxd';
denoised_cellset = convertStringsToChars(Info_Table.Denoised_CNMFe_Name);
denoised_cellset_ed = 'denoised_cellset_ed.isxd';
cnmfe_csvfile = 'cnmfe-cellset.csv';

%% 2.2: run preprocess

if ~exist(pp_file, 'file') && ~exist(cnmfe_file, 'file')
    isx.preprocess(...
        raw_movie_file, pp_file, 'temporal_downsample', Info_Table.temp_ds, ...
        'spatial_downsample', Info_Table.spat_ds, ...
        'fix_defective_pixels', true, ...
        'trim_early_frames', true);
end

%% Step 3: Spatial bandpass filtering: to increase signal contrast
%% 3.1: apply bandpass filter
% apply bandpass filter
if ~exist(bp_file, 'file') && ~exist(cnmfe_file, 'file')
    isx.spatial_filter(...
        pp_file, bp_file, ...
        'low_cutoff', 0.005, ...
        'high_cutoff', 0.5, ...
        'retain_mean', false);
end

%% Step 4: Motion correction: to stabilize the movie
%% 4.1 create reference image using temporal movie projection.

% project movie to get an mean image
if ~exist(mc_ref, 'file') && ~exist(cnmfe_file, 'file')
    isx.project_movie(bp_file, mc_ref, 'stat_type', 'mean');
end

%% 4.2: draw polygon as ROI for motion correction

% This defines the ROI used for motion correction, which is around the
% center of where cells usually show up (no need to change unless cells are REALLY off)

% if you want to change the coordinates, go to IDPS and process until this
% step, check the XY coordinates (shown around bottom of the IDPS
% interface) and enter them here
polygon_coordinates = [60,63;60,160;150,160;160,60];

%% 4.3 make motion corrected filename and two csv filenames 

% to store xy shift and corrected image rectangle reference
output_translation_csv = 'motion-corr-translation.csv';
output_crop_rect_csv = 'motion-corr-region.csv';

%% 4.4 run motion correction

% really no need to change any parameters here
if ~exist(mc_file, 'file') && ~exist(cnmfe_file, 'file')
    isx.motion_correct(...
        bp_file, mc_file, 'max_translation', 20, ...
        'low_bandpass_cutoff', 0, 'high_bandpass_cutoff', 0, ...
        'roi', polygon_coordinates, ...
        'reference_file_name', 'mc_ref.isxd', ...
        'global_registration_weight', 1, ...
        'output_translation_files', output_translation_csv, ...
        'output_crop_rect_file', output_crop_rect_csv);
end

%% Step 5: Run deltaF/F0 algorithm to get dff movie file for visualization & max image projection
%% 5.1 calculate dF/F0 & project the movie

% generate dF/F0 file
if ~exist(dff_file, 'file') && ~exist(cnmfe_file, 'file')
    isx.dff(mc_file, dff_file, 'f0_type', 'mean');
end

% maximum projection of the dff movie
if ~exist(maxdff_file, 'file') && ~exist(cnmfe_file, 'file')
    isx.project_movie(dff_file, maxdff_file, 'stat_type', 'max');
end


%% Step 6: Use CNMFe algorithm to segment cells 

%% 6.1 set parameters: refer to IDPS online documentation for "parameter setting tips"

% those parameters don't have to be changed for each animal, so keep them
% here
bg_spatial_subsampling=2;
ring_size_factor=1.4;
gaussian_kernel_size=0; % this is set to 0 = auto in isx.cnmfe
closing_kernel_size=0; % this is set to 0 = auto in isx.cnmfe
merge_threshold=0.7;
processing_mode=2; % processing each recording in parallel mode
num_threads=4;
patch_size=80;
patch_overlap=20;
output_unit_type=1;

%% 6.2: run cnmfe

if ~exist(cnmfe_file, 'file')
    tic
    isx.cnmfe(...
        mc_file, cnmfe_file, project_folder, 'cell_diameter', Info_Table.cell_diameter, ...
       'min_corr', Info_Table.min_corr, 'min_pnr', Info_Table.min_pnr, ...
       'bg_spatial_subsampling', bg_spatial_subsampling, 'ring_size_factor', ring_size_factor, ...
       'gaussian_kernel_size', gaussian_kernel_size, 'closing_kernel_size', closing_kernel_size, ...
       'merge_threshold', merge_threshold, 'processing_mode', processing_mode, ...
       'num_threads', num_threads, 'patch_size', patch_size, 'patch_overlap', patch_overlap, ...
       'output_unit_type', output_unit_type);
    toc
end

%% Step 7: Deconvolve cell traces from CNMFe for better estimation of neural calcium spikes 
%% 7.1: define parameters 

% those parameters don't have to be changed for each animal, keep here
accepted_only = false;
noise_method = 0;
first_order_ar = true;
lags = 5;
fudge_factor = 0.96;
% if you want another method of deconvolution you could check
% isx.deconvolve_cellset and decide, but this seems to work the best
deconvolution_method = 0; 

%% 7.2: run cellset deconvolution & event detection

if ~isfile(fullfile(Info_Table.CNMFe_Directory,Info_Table.Denoised_CNMFe_Name))
    delete(denoised_cellset_ed);
    isx.deconvolve_cellset(...
        cnmfe_file, denoised_cellset, denoised_cellset_ed, ...
        'accepted_only', accepted_only, 'spike_snr_threshold', Info_Table.spike_snr_threshold, ...
        'noise_range_min', Info_Table.noise_range_min, 'noise_range_max', Info_Table.noise_range_max, ...
        'noise_method', noise_method, 'first_order_ar', first_order_ar, ...
        'lags', lags, 'fudge_factor', fudge_factor, 'deconvolution_method', deconvolution_method);
end

%% Step 8: auto accept or reject cells

% set filter based on cell size, SNR, event rate or # components.

% # of components>1 if a single cell trace corresponds to several spatial patches, 
% which is unreasonable at least for a hippocampal pyr-cell recording.

% cells that have less than 5 events over 10min will be rejected
events_filters = {{'# Comps', '=', 1}, {'Event Rate', '>', 0.0083333}};

% run auto classification
if exist(denoised_cellset, 'file')
    isx.auto_accept_reject(denoised_cellset, denoised_cellset_ed, 'filters', events_filters);
end

%% Step 9: export cell calcium traces for CNMFe dataset

% run export function
if ~exist(cnmfe_csvfile, 'file')
    isx.export_cell_set_to_csv_tiff(...
        denoised_cellset, cnmfe_csvfile, 'cnmfe-cell.tiff', 'time_ref', 'start');
end

% list all tiff files
file_list = dir('*.tiff');
% get the file names
file_name = {file_list(:).name}';

% make subfolder dir and move files

if ~exist('cnmfe_tiff', 'dir')
    mkdir cnmfe_tiff
    % set subfolder dir and move files
    subfolder = fullfile(project_folder, 'cnmfe_tiff');
    movetiff = cellfun(@(each_cell) movefile(each_cell, subfolder), file_name);
end

%% Step 10: delete files that are generated in intermediate steps 
% don't need them and they take up quite amount of computer space

delete('prepro.isxd');
delete('mc_ref.isxd');
delete('bp_filter.isxd');
delete('motion-corr.isxd');
delete('motion-corr-region.csv');
delete('motion-corr-translation.csv');

%% Step 11: transfer denoised_cnmfe_cellset into CNMFe folder for LR analysis

% copy the denoised-cellset to another directory to get ready for LR
% analysis there
if exist(denoised_cellset, 'file')
    copyfile(Info_Table.Denoised_CNMFe_Name, Info_Table.CNMFe_Directory);
end

end