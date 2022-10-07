function output = importXYData(workbookFile, sheetName, dataLines)
%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [39, 15039];
end

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "A" + dataLines(1, 1) + ":P" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["Trialtime", "Recordingtime", "Xcenter", "Ycenter", "Xnose", "Ynose", "Xtail", "Ytail", "Area", "Areachange", "Elongation", "Direction", "Distancemoved", "Headdirectedtozone", "Headdirectedtozone2", "Result1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
output = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "A" + dataLines(idx, 1) + ":P" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    output = [output; tb]; %#ok<AGROW>
end

%% Convert to output type
output = table2array(output);
end