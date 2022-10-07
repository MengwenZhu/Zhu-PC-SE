function varargout=sysInfo(doDisp)

% % MIT License
% % 
% % Copyright (c) [2022] [Mengwen Zhu]
% % 
% % Permission is hereby granted, free of charge, to any person obtaining a copy
% % of this software and associated documentation files (the "Software"), to deal
% % in the Software without restriction, including without limitation the rights
% % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % copies of the Software, and to permit persons to whom the Software is
% % furnished to do so, subject to the following conditions:
% % 
% % The above copyright notice and this permission notice shall be included in all
% % copies or substantial portions of the Software.
% % 
% % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % SOFTWARE.

% sysInfo Retrieves information on system.
% sysInfo retrieves information on some performance-related specs of the
% system's hardware (CPU, GPU), including snapshots of the memory currently
% available. Additionally, it determines the Matlab version as well as the
% OS and sets a 'isWin' flag (because some Matlab functionality exists only
% on win).
% If called without input and output arguments, it will display the
% information in the command window.
% info = sysInfo  returns the information in struct array info.
% info = sysInfo(true) additionally displays the gathered information in
% the command window.

r.infoText = "--------- System Info ---------";

% CPU cores (be cautious, as this involves undocumented 'feature'):
try
  % - retrieve number of workers
  r.numWorker = feature('numcores');
  % - generate accompanying info text
  coreInfo = evalc('feature(''numcores'')');
  % note that the two lines below hinge on the exact output of function feature
  coreInfo = strsplit(string(coreInfo), '\n');
  % find "ans =" and dump everything including and after it
  tmpIx = find(coreInfo == "ans =");
  r.infoText = vertcat(r.infoText, coreInfo(1:tmpIx-1)');
catch
  r.numWorker = NaN;
  r.infoText = vertcat(r.infoText, "Number of CPU cores could not be determined");
end

% RAM and maximal size of array in Mb
r.ram = NaN;
r.ramAvailable = NaN;
r.maxArraySize = NaN;
if ispc
  [userview, sysview] = memory;
  % total
  r.ram = sysview.PhysicalMemory.Total/2^20;
  % currently available
  r.ramAvailable = sysview.PhysicalMemory.Available/2^20;
  % maximal possible array in Matlab (is a function of both swap space and
  % memory fragmentation)
  r.maxArraySize = userview.MaxPossibleArrayBytes/2^20;
elseif isunix
  % feature('MemStats') doesn't work on all Unix/Linux distros; use 'free'
  s = strsplit(string(evalc('!free -m')));
  tmpIx = find(s == "Mem:");
  % we assume that total and available memory are the values in position 1
  % and 6, respectively, after keyword "Mem:"
  r.ram = str2double(s(tmpIx+1));
  r.ramAvailable = str2double(s(tmpIx+6));
  % max array size cannot be determined
elseif ismac
  % try feature('MemStats') on MAC
  try
    s = strsplit(string(evalc('feature(''MemStats'')')));
    tmpIx = find(s == "Free:", 1);
    tmpIx2 = find(s == "Total:", 1);
    % as in the Unix case, assume a fixed output format
    r.ram = str2double(s(tmpIx2+1));
    r.ramAvailable = str2double(s(tmpIx+1));
    % second call to get the output value, the max possible array
    r.maxArraySize = feature('MemStats');
  catch
    r.infoText = vertcat(r.infoText, ...
      "Amount of RAM could not be determined:");
  end
else 
  % this should never be the case, though
  error("OS is neither of Win, Linux or MAC")
end
r.infoText = vertcat(r.infoText, "RAM (Megabytes)");
r.infoText = vertcat(r.infoText, "   Total: " + int2str(r.ram));
r.infoText = vertcat(r.infoText, "   Available: " + int2str(r.ramAvailable));
r.infoText = vertcat(r.infoText, "Maximal size of array: " + int2str(r.maxArraySize));

% OS
r.os = computer;
r.isWin = ispc;
r.infoText = vertcat(r.infoText, ...
  "Operating system: " + r.os);

% Matlab
r.MatlabVersion = version;
r.infoText = vertcat(r.infoText, ...
  "Matlab version: " + r.MatlabVersion);

% before checking GPUs, we need to know whether the parallel computing
% toolbox is installed
product_info = ver('parallel');
if ~isempty(product_info) && strcmpi(product_info.Name,'parallel computing toolbox')
    r.hasParallelToolbox = true;
else
    r.hasParallelToolbox = false;
end

% GPU
if r.hasParallelToolbox
    if gpuDeviceCount >= 1
        try
            gpuDev = gpuDevice;
            r.hasGPU = true;
            r.gpuName = gpuDev.Name;
            r.gpuMem = gpuDev.TotalMemory/2^20;
            r.gpuMemAvailable = gpuDev.AvailableMemory/2^20;
        catch ME
            % !be a little more specific here
            r.hasGPU = false;
            r.gpuName = "None (you may have to update the driver)";
            r.gpuMem = 0;
            r.gpuMemAvailable = 0;
        end
        r.infoText = vertcat(r.infoText, "GPUs");
        r.infoText = vertcat(r.infoText, ...
            "   Number of available GPUs: " + string(int2str(gpuDeviceCount)));
        r.infoText = vertcat(r.infoText, ...
            "   Currently selected GPU: " + r.gpuName);
        r.infoText = vertcat(r.infoText, ...
            "   GPU memory (Megabytes)");
        r.infoText = vertcat(r.infoText, ...
            "      Total: " + int2str(r.gpuMem));
        r.infoText = vertcat(r.infoText, ...
            "      Available: " + int2str(round(r.gpuMemAvailable)));
    else
        r.infoText = vertcat(r.infoText, ...
            "No Matlab-compatible GPU detected");
    end
else
    r.gpuName = "Cannot be determined";
    r.gpuMem = 0;
    r.gpuMemAvailable = 0;
    r.infoText = vertcat(r.infoText, ...
        "Matlab Parallel Computing Toolbox not available");
end

if nargout == 0 || (nargin > 0 && doDisp>0)
  disp(char(r.infoText))
end
if nargout > 0
  varargout{1} = r;
end