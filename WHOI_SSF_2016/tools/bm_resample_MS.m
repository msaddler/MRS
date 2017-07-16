function bm_resample_MS()
% MS 2016.06.23
% Blue whale data downsampling script

% bm14: fs0 = 500kHz
% bm15: fs0 = 120kHz
% bm16: fs0 = 120kHz

% ds: [5,10,10]
% 500,000 -> 100,000 -> 10,000 -> 1000

% fs1 = 120000/4;
% fs2 = 120000/4/3;
fs3 = 500000/5/10/10; % <-- bm14 is ds [5,10,10], bm15-16 is ds [4 3 10]

%if isempty(baseDir)
    baseDir = 'G:\bm14\'; % Base directory of all deployments in a year
    saveDir = 'G:\bm_ds\'; % Output directory for all downsampled data
%end

% Get all of the deployment subfolders in baseDir
depDir = dir(fullfile(baseDir));
depIdx = [depDir(:).isdir];
depNames = {depDir(depIdx).name}';
depNames = depNames(~ismember(depNames,{'.','..'}));

% Loop through all deployments in baseDir
for itr0=1:length(depNames)
    depFold = [baseDir,char(depNames(itr0)),'\']; % Deployment folder
    wavFiles = dir(fullfile(depFold,'*.wav'));
    
    % Set up deployment folder in output directory
    saveDepDir = [saveDir,char(depNames(itr0)),'\'];
    if ~exist(saveDepDir,'dir')
        mkdir(saveDepDir);
    end
    
    fn = []; % Get all .wav file names in the deployment folder
    for itr1 = 1:length(wavFiles)
        fn{1,itr1} = [depFold,wavFiles(itr1).name];
    end
    
    % Load all .wav files, resample, and save
    for itr1 = 1:length(fn)
        y = []; fs0 = [];
        disp(['Processing file: ',fn{itr1}])
        [y,fs0] = audioread_singlechannelds(fn{itr1});
        disp(['fs0 = ',num2str(fs0), 'Hz'])
        
        % Multi-stage resampling
        %y1 = resample(y,fs1,fs0); % taken care of in audioread_singlechannelds
        %y2 = resample(y,fs2,fs1); % likwise
        y3 = resample(y,fs3,fs0);
        
        % Save file with latest sample rate
        disp(['Saving downsampled file: ',num2str(fs3),' Hz'])
        saveName = [saveDepDir,char(wavFiles(itr1).name)];
        audiowrite(saveName,y3,fs3)
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified audioread to downsample in fragments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,Fs] = audioread_singlechannelds(filename, range, datatype)
% MS 2016.06.24
%
% Modified version of audioread
% - designed to read high sample rate recordings in small chunks
% - returns recording downsampled by 20
% - 2-stage downsampling [4, 5]

chunk = 6000000; % Size of fragment to read at once

%AUDIOREAD Read audio files
%   [Y, FS]=AUDIOREAD(FILENAME) reads an audio file specified by the
%   string FILE, returning the sampled data in Y and the sample rate
%   FS, in Hertz.
%
%   [Y, FS]=AUDIOREAD(FILENAME, [START END]) returns only samples START
%   through END from each channel in the file.
%
%   [Y, FS]=AUDIOREAD(FILENAME, DATATYPE) specifies the data type format of
%   Y used to represent samples read from the file.
%   If DATATYPE='double', Y contains double-precision normalized samples.
%   If DATATYPE='native', Y contains samples in the native data type
%   found in the file.  Interpretation of DATATYPE is case-insensitive and
%   partial matching is supported.
%   If omitted, DATATYPE='double'.
%
%   [Y, FS]=AUDIOREAD(FILENAME, [START END], DATATYPE);
%
%   Output Data Ranges
%   Y is returned as an m-by-n matrix, where m is the number of audio
%   samples read and n is the number of audio channels in the file.
%
%   If you do not specify DATATYPE, or dataType is 'double',
%   then Y is of type double, and matrix elements are normalized values
%   between -1.0 and 1.0.
%
%   If DATATYPE is 'native', then Y may be one of several MATLAB
%   data types, depending on the file format and the BitsPerSample
%   of the input file:
%
%    File Format      BitsPerSample  Data Type of Y     Data Range of Y
%    ----------------------------------------------------------------------
%    WAVE (.wav)            8           uint8             0 <= Y <= 255
%                          16           int16        -32768 <= Y <= 32767
%                          24           int32         -2^32 <= Y <= 2^32-1
%                          32           int32         -2^32 <= Y <= 2^32-1
%                          32           single         -1.0 <= Y <= +1.0
%    ----------------------------------------------------------------------
%    FLAC (.flac)           8           uint8             0 <= Y <= 255
%                          16           int16        -32768 <= Y <= 32767
%                          24           int32         -2^32 <= Y <= 2^32-1
%    ----------------------------------------------------------------------
%    MP3 (.mp3)            N/A          single         -1.0 <= Y <= +1.0
%    MPEG-4(.m4a,.mp4)
%    OGG (.ogg)
%    ----------------------------------------------------------------------
%
%   Call audioinfo to learn the BitsPerSample of the file.
%
%   Note that where Y is single or double and the BitsPerSample is
%   32 or 64, values in Y might exceed +1.0 or -1.0.
%
%   See also AUDIOINFO, AUDIOWRITE

%   Copyright 2012-2013 The MathWorks, Inc.

% Parse input arguments:
narginchk(1, 3);

if nargin < 2
    range = [1 inf];
    datatype = 'double';
elseif nargin < 3 && ischar(range)
    datatype = range;
    range = [1 inf];
elseif nargin < 3
    datatype = 'double';
end


% Expand the path, using the matlab path if necessary
filename = audiovideo.internal.absolutePathForReading(...
    filename, ...
    'MATLAB:audiovideo:audioread:fileNotFound', ...
    'MATLAB:audiovideo:audioread:filePermissionDenied');

import multimedia.internal.audio.file.PluginManager;

try
    readPlugin = PluginManager.getInstance.getPluginForRead(filename);
    
    options.Filename = filename;
    
    % Create Channel object
    channel = asyncio.Channel( ...
        readPlugin,...
        PluginManager.getInstance.MLConverter, ...
        options, [0, 0]);
    
    channel.InputStream.addFilter( ...
        PluginManager.getInstance.TransformFilter, ...
        []);
    
    
    % Validate the datatype is correctly formed
    
    pieces = 1:chunk:channel.TotalSamples;
    if pieces(end) ~= channel.TotalSamples
        pieces = [pieces,channel.TotalSamples];
    end
    y = [];
    for i = 1:length(pieces)-1
        range = [pieces(i), pieces(i+1)];
        disp(['Progress: ', num2str(i) ,'/', num2str(length(pieces)-1)])
        
        datatype = validateDataType(datatype, channel);
        
        [startSample, samplesToRead] = validateRange(range, channel.TotalSamples);
        
        options.StartSample = startSample;
        options.FrameSize = double(multimedia.internal.audio.file.FrameSize.Optimal);
        options.FilterTransformType = 'DeinterleaveTranspose';
        options.FilterOutputDataType = datatype;
        
        channel.open(options);
        c = onCleanup(@()channel.close()); % close when going out of scope
        
        % Read section of samples and keep only channel 1
        y_tmp = channel.InputStream.read(samplesToRead);
        y_tmp = y_tmp(:,1);
        
        % Two-stage downsampling (ds [4, 3])
        fs0 = double(channel.SampleRate);
        fs1 = fs0/5;
        fs2 = fs0/5/10;
        y_tmp = resample(y_tmp,fs1,fs0);
        y_tmp = resample(y_tmp,fs2,fs1);
        y = [y;y_tmp];
        
    end
    Fs = fs2;
    
    % Generate a warning if the number of samples read is lesser than the
    % number requested.
    actNumSamplesRead = size(y, 1);
    if actNumSamplesRead < samplesToRead
        stopSample = startSample + samplesToRead - 1;
        if (startSample+1 ~= 1) || (stopSample+1 ~= channel.TotalSamples)
            warning(message('MATLAB:audiovideo:audioread:incompleteRead', actNumSamplesRead));
        end
    end
    
catch exception
    exception = PluginManager.convertPluginException(exception, ...
        'MATLAB:audiovideo:audioread');
    
    throw(exception);
end

end

function datatype = validateDataType(datatype, channel)

datatype = validatestring(datatype, {'double','native'},'audioread','datatype');

if strcmp(datatype,'native')
    if ismember('BitsPerSample',properties(channel))
        % Channel has a 'BitsPerSample' property and is most likely
        % uncompressed or lossless.
        % Set the 'native' data type to the underlying channel's
        % datatype.
        datatype = channel.DataType;
    else
        % Channel is most likely compressed. 'native' datatype
        % should be single
        datatype = 'single';
    end
else
    datatype = 'double';
end

end

function [startSample, samplesToRead] = validateRange(range, totalSamples)

validateattributes( ...
    range,{'numeric'}, ...
    {'positive','nonempty','nonnan','ncols',2,'nrows',1}',...
    'audioread','range',2);

% replace any Inf values with total samples
range(range == Inf) = totalSamples;

if any(range > totalSamples)
    error(message('MATLAB:audiovideo:audioread:endoffile', totalSamples));
end


range = range - 1; % sample ranges are zero based
range = max(range, 0); % bound the range by zero

% Validate that all values are integers
validateattributes( ...
    range,{'numeric'}, ...
    {'integer'},...
    'audioread','range',2);

% Validate that start of range is less than the end of range
if (range(1) > range(2))
    error(message('MATLAB:audiovideo:audioread:invalidrange'));
end

startSample = range(1);
samplesToRead = (range(2) - startSample) + 1;

end