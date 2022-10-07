function filePaths = getNoldusRawDataPathFromMetadata(metadata, dirStub)

formatSpec = "Locating Noldus raw data for %s, %s, %s...";
disp(compose(formatSpec, metadata.animalName, metadata.exptDate1))

% append genotype info to dir stub
dirStub = fullfile(dirStub, getDirStubPartFromGenotype(metadata));
% construct animal-specific subdirectory, first part (genotype info)
subdir = getSubdirPartFromGenotype(metadata);
% append drug info to animal-specific subdirectory
% subdir = [subdir, '_', getSubdirPartFromDrug(metadata)];
% add subdirectory in which all CNMFe files are collected - this completes
% the path to the animal's raw CNMFe data files
subdir = fullfile(dirStub, subdir, 'Noldus Files');

if strcmp(metadata.exptParadigm, '4h_memory_test') == 1
    % look within Noldus folder for file corresponding to this date
    dateToFind = strrep(metadata.exptDate1,'-','_'); % replace the dashes (-) with underscores (_)
    
    % since the Noldus files have a single digit in the months place (but not
    % the Inscopix!), remove the zero digit from the months place in dateToFind
    underscorePos = strfind(dateToFind,'_'); % position of the underscore, should be after the date
    removeDigit = false(size(underscorePos)); % create logical array to use to remove these extra zeros
    for ii=1:size(underscorePos,2)
        if strcmp(dateToFind(underscorePos(ii)+1),'0')
           removeDigit(ii) = true;
        end
    end
    charsToDelete = underscorePos(removeDigit)+1; % remove the digit AFTER the underscore
    dateToFind(charsToDelete) = '';

    D = dir(subdir);
    listOfFiles = {D.name};

    % am file
    am = contains(listOfFiles,dateToFind) & contains(listOfFiles,'am','IgnoreCase',true);
    amFile = listOfFiles{am};
    % pm file
    pm = contains(listOfFiles,dateToFind) & contains(listOfFiles,'pm','IgnoreCase',true);
    pmFile = listOfFiles{pm};
    % get directory to these noldus files
    filePaths.AM = fullfile(subdir, amFile);
    filePaths.PM = fullfile(subdir, pmFile);
    
elseif strcmp(metadata.exptParadigm, '24h_memory_test') == 1
    dateToFind1 = strrep(metadata.exptDate1,'-','_'); 
    dateToFind2 = strrep(metadata.exptDate2,'-','_');
    % for date 1
    underscorePos = strfind(dateToFind1,'_'); 
    removeDigit = false(size(underscorePos)); 
    for ii=1:size(underscorePos,2)
        if strcmp(dateToFind1(underscorePos(ii)+1),'0')
           removeDigit(ii) = true;
        end
    end
    charsToDelete = underscorePos(removeDigit)+1; 
    dateToFind1(charsToDelete) = '';
    % for date 2
    underscorePos = strfind(dateToFind2,'_'); % position of the underscore, should be after the date
    removeDigit = false(size(underscorePos)); % create logical array to use to remove these extra zeros
    for ii=1:size(underscorePos,2)
        if strcmp(dateToFind2(underscorePos(ii)+1),'0')
           removeDigit(ii) = true;
        end
    end
    charsToDelete = underscorePos(removeDigit)+1; % remove the digit AFTER the underscore
    dateToFind2(charsToDelete) = '';

    D = dir(subdir);
    listOfFiles = {D.name};

    % day1 file
    day1 = contains(listOfFiles,dateToFind1) & contains(listOfFiles,'day1','IgnoreCase',true);
    day1File = listOfFiles{day1};
    % day2 file
    day2 = contains(listOfFiles,dateToFind2) & contains(listOfFiles,'day2','IgnoreCase',true);
    day2File = listOfFiles{day2}; 
    % get directories to these noldus files
    filePaths.day1 = fullfile(subdir, day1File);
    filePaths.day2 = fullfile(subdir, day2File);
    
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
