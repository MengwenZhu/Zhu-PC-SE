function [CNMFefilePaths, filePaths, fileNames] = getCaImagingRawRecordingPathAndNameFromMetadata(metadata, dirStub)
% getCaImagingRawRecordingPathAndNameFromMetadata - construct full paths to
% a pair of am/pm or day1/day2 calcium imaging raw recording, full paths to 
% the corresponding CNMFe folder, and full name of the recordings
% 
% [CNMFefilePaths, filePaths, fileNames] = getCaImagingRawRecordingPathFromMetadata(metadata, dirStub)
% constructs full path to a pair of am/pm or day1/day2 raw calcium imaging recordings. 
% Input arg metadata is a struct with fields animalName,
% exptDate, drug, dose, dataset, genotype, and exptParadigm: it is a single row of the
% high-level metadata file used for collecting all experiments, converted
% to a struct. Input arg dirStub is the base path to where calcium imaging recordings are organized 
% (ds.CNMFe_Inscopix_rawdataPath).
% Output arg CNMFefilePaths is a structure with directory of a folder named "CNMFe Files"
% that is used to collect denoised_cnmfe_cellsets for that animal
% filePaths is a structure with a pair of paths to the pair of recordings
% fileNames is a structure with a pair of recording names

formatSpec = "Locating raw calcium imaging recording for %s, %s, %s...";
disp(compose(formatSpec, metadata.animalName, metadata.exptDate1))

% append genotype info to dir stub
dirStub = fullfile(dirStub, getDirStubPartFromGenotype(metadata));
% construct animal-specific subdirectory, first part (genotype info)
subdir = getSubdirPartFromGenotype(metadata);
% add subdirectory in which all folders containing the raw recordings
% locate
subdir = fullfile(dirStub, subdir);


% now construct patterns of expected calcium imaging recording directories & recording names 
% you will get a pair of file names and directories (AM&PM or day1&day2)
if strcmp(convertCharsToStrings(metadata.exptParadigm), '4h_memory_test') == 1
    filePaths = struct('AM', string(missing), 'PM', string(missing));
    tod = string(fieldnames(filePaths));
    for ix = 1:numel(tod)
        fileNamePattern = metadata.exptDate1 + " " + tod(ix) + "*";
        ff = dir(fullfile(subdir, fileNamePattern));
        if length(ff) == 1
            filePaths.(tod(ix)) = string(fullfile(subdir, ff.name));
        elseif isempty(ff)
            warning("No " + tod(ix) + " file found");
        else
            warning("Several " + tod(ix) + " files found");     %TODO: list them
        end
    end
elseif strcmp(convertCharsToStrings(metadata.exptParadigm), '24h_memory_test') == 1
     filePaths = struct('day1', string(missing), 'day2', string(missing));
     filePaths.day1 = metadata.exptDate1 + " " + 'day1' + "*";
     ff = dir(fullfile(subdir, filePaths.day1));
     filePaths.day1 = string(fullfile(subdir,ff.name));
     filePaths.day2 = metadata.exptDate2 + " " + 'day2' + "*"; 
     ff = dir(fullfile(subdir, filePaths.day2));
     filePaths.day2 = string(fullfile(subdir,ff.name));
end


% Now generate recording names (raw calcium triggered recording) 
if strcmp(convertCharsToStrings(metadata.exptParadigm), '4h_memory_test') == 1
    fileNames = struct('AM', string(missing), 'PM', string(missing));
    tod = string(fieldnames(fileNames));
    for ix = 1:numel(tod)
        ff = dir(filePaths.(tod(ix)));
        for kk=1:height(ff)          
            if contains(ff(kk).name, 'video_trig')
                fileNames.(tod(ix)) = ff(kk).name;
            end
        end
    end
elseif strcmp(convertCharsToStrings(metadata.exptParadigm), '24h_memory_test') == 1
    fileNames = struct('day1', string(missing), 'day2', string(missing));
    tod = string(fieldnames(fileNames));
    for ix = 1:numel(tod)
        ff = dir(filePaths.(tod(ix)));
        for kk=1:height(ff)  
            if contains(ff(kk).name, 'video_trig')
                fileNames.(tod(ix)) = ff(kk).name;
            end
        end
    end    
else
    error('currently, only expecting 4h and 24h memory tests')
end
    

% create CNMFe file Path
CNMFefilePaths = struct();
CNMFefilePaths.path = fullfile(convertCharsToStrings(subdir), 'CNMFe Files');


% ---------------------- LOCAL FUNCTIONS ----------------------------------

function dirStub = getDirStubPartFromGenotype(metadata)
% getDirStubPartFromGenotype - amend dirStub according to genotype info
if contains(metadata.animalName, 'G2')
    % define subdirectory indicating genotype type
    dirStub = 'GAD Mice';
elseif contains(metadata.animalName, 'C57')
    dirStub = 'C57 Mice';
elseif contains(metadata.animalName, 'CCK')
    dirStub = 'CCK Mice';
else
    error('currently, expecting only GAD, C57, and CCK mice')
end

function subdir = getSubdirPartFromGenotype(metadata)
% getSubdirPartFromGenotype - construct first part of animal's expected
% subdirectory name based on genotype
if contains(metadata.animalName, 'G2')
    subdir = ['GAD-', strrep(strrep(metadata.animalName, 'G2_', ''), '_', '-'), ' '];
    % 'KO' (knockout) translates to 'Cre+', 'WT' to 'Cre-'
    if contains(metadata.genotype, 'KO')
        subdir = [subdir, 'Cre+'];
    elseif contains(metadata.genotype, 'WT')
        subdir = [subdir, 'Cre-'];
    else
        error('missing information on genotype');
    end
elseif contains(metadata.animalName, 'C57')
    subdir = [metadata.animalName];
elseif contains(metadata.animalName, 'CCK')
    subdir = [metadata.animalName];
    if contains(metadata.genotype, 'KO')
        subdir = [subdir, ' ', 'Cre+'];
    elseif contains(metadata.genotype, 'WT')
        subdir = [subdir, ' ', 'Cre-'];
    elseif contains(metadata.genotype, 'double-KO')
        subdir = [subdir, ' ', 'Cre++'];
    else
        error('missing information on genotype');
    end
else
    error('currently, expecting only GAD, C57, and CCK mice')
end
