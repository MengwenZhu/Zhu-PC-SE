function filePaths = getInscopixRawDataPathFromMetadata(metadata, dirStub)

formatSpec = "Locating Inscopix raw data for %s, %s, %s...";
disp(compose(formatSpec, metadata.animalName, metadata.exptDate1))

% append genotype info to dir stub
dirStub = fullfile(dirStub, getDirStubPartFromGenotype(metadata));
% construct animal-specific subdirectory, first part (genotype info)
subdir = getSubdirPartFromGenotype(metadata);
% add subdirectory in which all CNMFe files are collected - this completes
% the path to the animal's raw CNMFe data files
subdir = fullfile(dirStub, subdir, 'CNMFe Files', 'LR export events');

% now construct patterns of expected file names
if strcmp(convertCharsToStrings(metadata.exptParadigm), '4h_memory_test') == 1
    filePaths = struct('AM', string(missing), 'PM', string(missing));
    tod = string(fieldnames(filePaths));
    for ix = 1:numel(tod)
        fileNamePattern = metadata.exptDate1 + " " + tod(ix) + "*.csv";
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
    filePaths.day1 = metadata.exptDate1 + " " + 'day1' + "*.csv";
    ff = dir(fullfile(subdir, filePaths.day1));
    filePaths.day1 = string(fullfile(subdir,ff.name));
    filePaths.day2 = metadata.exptDate2 + " " + 'day2' + "*.csv";
    ff = dir(fullfile(subdir, filePaths.day2));
    filePaths.day2 = string(fullfile(subdir,ff.name));
end


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
