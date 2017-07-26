function dtagmark()
% MS 2016.06.17
%
% Custom graphic user interface designed for auditing low-frequency blue
% whale calls, especially for aligning audio recordings with accelerometer
% readings.
%
% KNOWN BUG: saving audit log to .xlsx file relies on "xlswrite" MATLAB
% function, which causes error in MAC OSX.
%

clear; close all;

dataDir = '/Users/mark/Documents/WHOI_SSF_2016/Mark/data_bm/';

%%% PARAMETERS
TW = 60; % Time window (s)
TW_CAP = 180; % Cap the allowed size of TW to prevent crash (s)
pt = 256; % FFT point size for audio spectrogram
pt_a = 256; % FFT point size for accelerometer spectrogram
overlap = 0.90; % FFT fraction of overlap for audio spectrogram
overlap_a = 0.90; % FFT fraction of overlap for accelerometer spectrogram

CLIM = [-90, 15]; % Audio colorbar limits (dB)
FLIM = [0 500]; % Audio spectrogram frequency limits (Hz)

bandpass_a = [0,0]; % Accelerometer high-pass filter cutoff frequency (Hz)

%%% INITIALIZE VARIABLES
recDir = []; % Directory with audio (.wav) and metadata (.xml) files

fnIdx = 0; % Index of current audio file
fn = []; % Audio file names
fn_a = []; % Accelerometer file name

ts = 0; % Current window start time (MATLAB datenum, days)
T = []; % Cell array of time vectors for audio data (MATLAB datenums)
T_a = []; % Cell array of time vectors for accel. data (MATLAB datenums)

xmlTimes = 0; % Audio file start times (MATLAB datenums)
fs = 0; % Audio sample rate (Hz)
y = []; % Audio signals
Afs = 0; % Accelerometer sample rate (Hz)
A = []; % Accelerometer signals

audiospec = gobjects(1,1);
accelspec = gobjects(3,1);
waveFlag = 0; % Plot spectrogram (0) or waveform (1)

fn_log = []; % Audit log file name
log = {'DATENUM','DURATION','FILE','DATESTR','CALL','COMMENT'};
auditList = []; % Short form of log entries to display on GUI
commLength = 24; % Number of comment characters to display in log entries


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Graphic User Interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Name','DTAGmark','Position',[50,100,1200,800]);
pan_c = uipanel('Parent',fig,'Position',[.01 .01 .18 .98],'Title','Controls');
pan_p = uipanel('Parent',fig,'Position',[.20 .01 .79 .98],'Title','Plot');

% Control Panel Fields
c_bck = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .01 .48 .03],...
    'String','<','Callback',@callback_pushbutton);
c_fwd = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.50 .01 .48 .03],...
    'String','>','Callback',@callback_pushbutton);

c_loadaudio = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .96 .96 .03],...
    'String','Load Audio Directory','Callback',@callback_pushbutton);

uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .91 .58 .03],...
    'String','Deployment Time:');
c_ts = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .91 .38 .03],...
    'String',datestr(ts-xmlTimes(1),'HH:MM:SS'),'Callback',@callback_edit);
uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .86 .58 .03],...
    'String','Window Size (s):');
c_TW = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .86 .38 .03],...
    'String',num2str(TW),'Callback',@callback_edit);
uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .81 .58 .03],...
    'String','Audio FFT size:');
c_pt = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .81 .38 .03],...
    'String',num2str(pt),'Callback',@callback_edit);
uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .76 .58 .03],...
    'String','Accel FFT size:');
c_pt_a = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .76 .38 .03],...
    'String',num2str(pt_a),'Callback',@callback_edit);
uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .71 .58 .03],...
    'String','Audio % Overlap:');
c_overlap = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .71 .38 .03],...
    'String',num2str(overlap),'Callback',@callback_edit);
uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .66 .58 .03],...
    'String','Accel % Overlap:');
c_overlap_a = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .66 .38 .03],...
    'String',num2str(overlap_a),'Callback',@callback_edit);

c_loadaudit = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .60 .96 .03],...
    'String','Load Audit Log','Callback',@callback_pushbutton);
c_auditlist = uicontrol('Parent',pan_c,'Style','listbox',...
    'Units','Normalized','Position',[.02 .33 .96 .26],...
    'String',auditList,'Callback',@callback_listbox);

c_loadacc = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .27 .96 .03],...
    'String','Load PRH File','Callback',@callback_pushbutton);

c_plotselect = uibuttongroup('Parent',pan_c,...
    'Position',[.02 .20 .96 .05],'SelectionChangedFcn',@plotselect);
c_spec = uicontrol(c_plotselect,'Style','radiobutton',...
    'Units','Normalized','Position',[.05 .01 .44 .98],...
    'String','Spectrogram');
c_wave = uicontrol(c_plotselect,'Style','radiobutton',...
    'Units','Normalized','Position',[.55 .01 .44 .98],...
    'String','Waveform');

uicontrol('Parent',pan_c,'Style','Text',...
    'Units','Normalized','Position',[.02 .15 .58 .03],...
    'String','Accel band-pass filter (Hz):');
c_bandpass_a1 = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.60 .15 .19 .03],...
    'String',num2str(bandpass_a(1)),'Callback',@callback_edit);
c_bandpass_a2 = uicontrol('Parent',pan_c,'Style','edit',...
    'Units','Normalized','Position',[.79 .15 .19 .03],...
    'String',num2str(bandpass_a(2)),'Callback',@callback_edit);

% Plot Panel Fields
AXs = axes('Parent',pan_p,'position',[.05 .62 .90 .35]); % Audio axes
AXa = gobjects(3,1); % Accelerometer axes (3 channels)
AXa(1) = axes('Parent',pan_p,'position',[.05 .40 .90 .17]);
AXa(2) = axes('Parent',pan_p,'position',[.05 .22 .90 .17]);
AXa(3) = axes('Parent',pan_p,'position',[.05 .04 .90 .17]);
linkaxes([AXs;AXa],'x')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function callback_pushbutton(source,~)
        switch source
            case c_bck
                ts = ts - TW/(60*60*24); % Step back one window length
                plot_audiospectrogram()
            case c_fwd
                ts = ts + TW/(60*60*24); % Step forward one window length
                plot_audiospectrogram()
            case c_loadaudio
                % Request user to select directory with .wav and .xml files
                recDir = uigetdir(dataDir,...
                    'Select directory with .wav and .xml files');
                if recDir == 0
                    disp('Load operation canceled'); return
                end
                recDir = [recDir,'/'];
                % Get names of .wav and .xml files
                wavFiles = dir(fullfile(recDir,'*.wav'));
                xmlFiles = dir(fullfile(recDir,'*.xml'));
                % Check directory contains .xml and .wav files
                if isempty(wavFiles) || isempty(xmlFiles)
                    warning('Invalid Directory: .wav and .xml files required')
                    return
                elseif length(wavFiles) ~= length(xmlFiles)
                    warning('Invalid Directory: # .xml ~= # .wav files')
                    return
                end
                % Reset xmlTimes vector upon loading new deployment data
                xmlTimes = zeros(length(wavFiles),1);
                % Read-in start times of each chip from .xml files
                for itr0 = 1:length(wavFiles)
                    tmpXML = [recDir,xmlFiles(itr0).name]; % .xml file name
                    xmlTimes(itr0,1) = extractxmlstart(tmpXML); % local fcn
                end
                % Reset filenames and signals array upon loading new
                % deployment data
                fn = cell(length(wavFiles),1);
                y = cell(length(wavFiles),1);
                fs = 0;
                T = []; % Initialize cell array for time vectors
                % Read-in audio signals from .wav files
                for itr0 = 1:length(wavFiles)
                    fn{itr0} = [recDir,wavFiles(itr0).name];
                    y{itr0} = [];
                    % Check .wav file to ensure manageable sample rate
                    [~,fs_tmp] = audioread(fn{itr0},1:2);
                    if fs_tmp > 2000 % Check sample rate is <2kHz
                        warning('fs too large to open'); return
                    elseif itr0>1 && fs_tmp~=fs
                        warning('Sample rates are not consistent'); return
                    end
                    disp(['Reading: ',fn{itr0}, ' (file ',num2str(itr0),...
                        ' of ',num2str(length(wavFiles)),')'])
                    [y{itr0},fs] = audioread(fn{itr0});
                    % Set up time vectors for audio files (MATLAB datenums)
                    T{itr0} = 1:(length(y{itr0}));
                    T{itr0} = (T{itr0}/fs)/(60*60*24) + xmlTimes(itr0);
                end
                ts = xmlTimes(1); % Default current time to initial time
                plot_audiospectrogram()
                
            case c_loadacc
                fn_a = []; A = []; Afs = 0;
                if isempty(recDir) % Check audio directory is loaded first
                    warning('Load audio directory first'); return
                end
                % Determine valid prh file name from recDir
                aTag = ['*',recDir(end-9:end-1),'prh.mat'];
                [FileName,PathName] = uigetfile([dataDir,aTag],...
                    'Select prh.mat file to load');
                fn_a = [PathName,FileName];
                if ~ischar(fn_a) % Check file was selected
                    disp('Load operation canceled');
                    fn_a = []; return
                end
                disp(['Loading accelerometer data: ',fn_a])
                tmpS = load(fn_a,'A','fs');        % Load acc data (A) and
                A_tmp = tmpS.A; Afs = tmpS.fs;     % acc sample rate (Afs)
                % (1) Set up time vectors for accelerometer data (T_a) to
                % align with the time vectors for audio data (T)
                % (2) Set up cell array for accelerometer data (A) to align
                % with the audio data (y)
                T_a = []; % Initialize cell array for acc time vectors
                for itr0 = 1:length(fn)
                    t_start = xmlTimes(itr0);
                    t_end = T{itr0}(end); % Align T_a with audio times
                    T_a{itr0} = t_start:(1/Afs)/(60*60*24):t_end;
                    
                    A_end = min(length(A_tmp),length(T_a{itr0}));
                    A{itr0} = A_tmp(1:A_end,:);
                    A_tmp = A_tmp((length(T_a{itr0})+1):end,:);
                end
                plot_audiospectrogram()
            case c_loadaudit
                % Request user to open log audit (.xlsx) file
                [FileName,PathName] = uigetfile([dataDir,'*.xlsx'],...
                    'Select audit log (Excel spreadsheet)');
                fn_log = [PathName,FileName];
                % If no log opened, cancel load operation
                if ~ischar(fn_log)
                    disp('Load operation canceled'); return
                else
                    % Otherwise, read-in the selected log file to cell array
                    disp(['Loading audit log file: ',fn_log])
                    [num,txt,raw] = xlsread(fn_log);
                    headers = txt(1,:);
                    if sum(strcmp(headers,log(1,:)))~=6
                        warning('IMPROPER LOG FORMAT'); return
                    end
                    log = raw; % log is stored in a cell array
                    if ~isempty(num)
                        % Display log entries in short form on in listbox
                        auditList = num(:,1);
                        auditList = datestr(auditList,'yyyy/mm/dd HH:MM:SS');
                        calltype = char(txt(:,5));
                        calltype = calltype(2:end,1:4);
                        space = char(zeros(size(calltype,1),2));
                        space(:,:) = ' ';
                        txt{1,6} = [txt{1,6},char(zeros(1,commLength))];
                        comm = char(txt(:,6));
                        comm = comm(2:end,1:commLength);
                        auditList = [auditList,space,calltype,space,comm];
                        set(c_auditlist,'String',auditList,'Value',size(auditList,1))
                    end
                end
        end
    end

    function callback_edit(source,~)
        switch source
            case c_ts
                val = datevec(get(source,'String'));
                val(1:3) = 0; % Get 'HH:MM:SS' format
                val = datenum(val); % Get MATLAB datenum format
                ts = xmlTimes(1)+val; % Convert deployment time to start time
                plot_audiospectrogram()
            case c_TW
                val = str2double(get(source,'String'));
                if isnan(val)
                    warning('Numerical value required'); return
                end
                TW = min(abs(val),TW_CAP);
                set(c_TW,'String',num2str(TW))
                plot_audiospectrogram()
            case c_pt
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    pt = abs(val);
                    plot_audiospectrogram()
                end
            case c_pt_a
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    pt_a = abs(val);
                    plot_audiospectrogram()
                end
            case c_overlap
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    overlap = abs(val);
                    plot_audiospectrogram()
                end
            case c_overlap_a
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    overlap_a = abs(val);
                    plot_accelspectrogram()
                end
            case c_bandpass_a1
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    bandpass_a(1) = abs(val);
                    plot_accelspectrogram()
                end
            case c_bandpass_a2
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    bandpass_a(2) = abs(val);
                    plot_accelspectrogram()
                end
        end
    end

    function plotselect(source,event)
        switch event.NewValue
            case c_spec
                waveFlag = 0; % Indicate to plot spectrograms
                plot_audiospectrogram()
            case c_wave
                waveFlag = 1; % Indicate to plot waveforms
                plot_audiospectrogram()
        end
    end

    function callback_listbox(source,event)
        switch source
            case c_auditlist
                entryIdx = get(source,'Value');
                if ~isempty(entryIdx)
                    entryStr = auditList(entryIdx,:);
                    ts = datenum(entryStr(1:19),'yyyy/mm/dd HH:MM:SS');
                    plot_audiospectrogram()
                end
        end
    end

    function buttondown_log(source,event)
        if event.Button ~= 3 % Right-click to make log entry
            disp('Right-click to begin log entry'); return
        end
        [t_tmp,~] = ginput(2);
        call_t = t_tmp(1);
        call_dur = t_tmp(2)-t_tmp(1);
        if call_dur <= 0
            warning('Invalid selection'); return
        end
        % Request user to input call type and a comment (optional)
        prompt = {'Call Type:','Comment:'};
        response = inputdlg(prompt,'Log Audit');
        call_type = response{1};
        call_comment = response{2};
        % If a call type is provided, generate log entry
        if ~isempty(call_type)
            row = size(log,1)+1;
            log{row,1} = num2str(call_t,12);
            log{row,2} = num2str(call_dur,7);
            log{row,3} = fn{fnIdx};
            log{row,4} = datestr(call_t, 'yyyy/mm/dd HH:MM:SS');
            log{row,5} = call_type;
            log{row,6} = call_comment;
            xlswrite(fn_log,log) % Write log to the loaded log file
            % Update the listbox with the newest log entry
            call_type = [call_type,'    '];
            call_comment = [call_comment,char(zeros(1,commLength))];
            newLine = [log{row,4},'  ',call_type(1:4)...
                ,'  ',call_comment(1:commLength)];
            auditList = [auditList;newLine];
            set(c_auditlist,'String',auditList)
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_audiospectrogram()
        % Use current window start time (ts) to determine current audio
        % file
        
        % Check that audio (.wav) files have been loaded
        if isempty(fn)
            return
        end
        
        % Check that window start time is in-bounds
        if ts < xmlTimes(1)
            ts = xmlTimes(1);
        end
        
        fnIdx = find(xmlTimes <= ts,1,'last'); % Get file name index
        tsi = 1 + floor((ts-xmlTimes(fnIdx))*fs*60*60*24); % Window start time index
        tei = tsi+(TW*fs); % Window end time index
        
        % If window start index is right at the end of a file, jump to the
        % next wavfile (unless already at the end of the recording)
        if abs(tsi-length(y{fnIdx}))<=fs && fnIdx ~= length(xmlTimes)
            disp('Moving to next audio file')
            fnIdx = fnIdx + 1;
            tsi = 1;
            tei = tsi+(TW*fs);
            ts = xmlTimes(fnIdx) + ((tsi-1)/fs)/(60*60*24);
        end
        
        % Check that window end time is in-bounds
        if tei > length(y{fnIdx})
            warning('End of audio file reached')
            tsi = length(y{fnIdx})-fs*TW; % Choose last in-bounds window
            tei = length(y{fnIdx});       % start index
            % Convert start index into corrected window start time
            ts = xmlTimes(fnIdx) + (tsi/fs)/(60*60*24);
        end
        
        % Update the displayed deployment time with the correct window
        % start time (adjusted if proposed window is out of bounds)
        set(c_ts,'String',datestr(ts-xmlTimes(1),'HH:MM:SS'))
        
        yIdx = tsi:tei; % Grab signal indexes within time window
        y_tmp = y{fnIdx};
        y_tmp = y_tmp(yIdx);
        
        if waveFlag % If waveFlag == 1, plot waveform instead
            plot_audiowaveform(y_tmp)
            plot_accelspectrogram()
            return
        end
        
        % Calculate spectrogram
        [s,f,t] = spectrogram(y_tmp,hamming(pt),floor(overlap*pt),pt,fs);
        s = 10*log10(abs(s));
        
        % Plot spectrogram and format axis
        axes(AXs)
        audiospec = imagesc(ts+(t/(60*60*24)),f,s,...
            'ButtonDownFcn',@buttondown_log);
        colormap('bone')
        set(AXs,'Ydir','Normal','ylim',FLIM,'clim',CLIM)
        datetick('x','HH:MM:SS','keeplimits','keepticks')
        % Add title with audio filename and time-into-file
        title(AXs,['Window Start: ',datestr(ts),'    (',fn{fnIdx},...
            '   -   ',datestr(ts-xmlTimes(fnIdx),'dd / HH:MM:SS'),')'],...
            'FontSize',8,'Interpreter','none')
        
        plot_accelspectrogram()
    end

    function plot_accelspectrogram()
        % Plot accelerometer data spectrogram (aligned with audio spec)
        if isempty(fn_a)
            return
        end
        
        % Determine accelerometer sample indexes
        tsA = 1 + floor((ts-xmlTimes(fnIdx))*Afs*60*60*24); % Window start time index
        teA = tsA+(TW*Afs); % Window end time index
        A_tmp = A{fnIdx}(tsA:teA,:);
        
        % If cutoff frequencies are provided, band-pass filter accelerometer data
        if 0<bandpass_a(1) && bandpass_a(1)<Afs
            if bandpass_a(2) <= bandpass_a(1)
                [b,a] = butter(5,bandpass_a(1)/(Afs/2),'high');
                A_tmp = filter(b,a,A_tmp);
            elseif bandpass_a(1)<bandpass_a(2) && bandpass_a(2)<Afs
                [b,a] = butter(5,bandpass_a/(Afs/2),'bandpass');
                A_tmp = filter(b,a,A_tmp);
            end
        end
        
        if waveFlag % If waveFlag == 1, plot waveform instead
            plot_accelwaveform(A_tmp)
            return
        end
        
        for i = 1:length(AXa)
            % Calculate accelerometer data spectrogram
            [s,f,t] = spectrogram(A_tmp(:,i),...
                hamming(pt_a),floor(overlap_a*pt_a),pt_a,Afs);
            s = 10*log10(abs(s)); % Convert to decibels
            
            % Plot spectrogram and format axis
            axes(AXa(i))
            accelspec(i) = imagesc(ts+(t/(60*60*24)),f,s);
            colormap(AXa(i),'bone')
            set(AXa(i),'YDir','Normal')%,'ylim',FLIM_a,'clim',CLIM_a)
            datetick('x','HH:MM:SS','keeplimits','keepticks')
            if i~= length(AXa)
                set(AXa(i),'xtick',[]) % Hide x-axis labels
            end
        end
    end

    function plot_audiowaveform(y_tmp)
        % Plot audio recording waveform (called if waveFlag == 1)
        % y_tmp contains full window of audio recording (1 channel)
        t = 0:(1/fs):(length(y_tmp)-1)/fs;
        t = ts + t/(60*60*24);
        
        plot(AXs,t,y_tmp)
        set(AXs,'xlim',[t(1),t(end)])
        datetick(AXs,'x','HH:MM:SS','keeplimits','keepticks')
        
        % Add title with audio filename and time-into-file
        title(AXs,['Window Start: ',datestr(ts),'    (',fn{fnIdx},...
            '   -   ',datestr(ts-xmlTimes(fnIdx),'dd / HH:MM:SS'),')'],...
            'FontSize',8,'Interpreter','none')
    end

    function plot_accelwaveform(A_tmp)
        % Plot accelerometer waveform (called if waveFlag == 1)
        % A_tmp contains full window of accelerometer readings (3 channels)
        t = 0:(1/Afs):(length(A_tmp)-1)/Afs;
        t = ts + t/(60*60*24);
        
        for i = 1:length(AXa)
            plot(AXa(i),t,A_tmp(:,1))
            set(AXa(i),'xlim',[t(1),t(end)])
            datetick(AXa(i),'x','HH:MM:SS','keeplimits','keepticks')
            if i~= length(AXa)
                set(AXa(i),'xtick',[]) % Hide x-axis labels
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function output = extractxmlstart(XML)
        % Function to extract the recording start and end times from the
        % metadata stored in .XML file format
        %
        % Input: .XML file name or XML string
        % Output: [start datenum, end datenum]
        %
        % Modified from xml2mat.m (Jonas Almeida, 20 Aug 2002, XML4MAT Tbox)
        %
        % Check if input is an XML file name or an XML string
        if XML(1)~='<'
            XML=strrep(file2str(XML),'''','''''');
            % Remove non-content lines if they exist
            XML=regexprep(XML,'<[?!].*?>','');
        end
        % Analyze XML line
        tag_ini=find(XML=='<');
        tag_end=find(XML=='>');
        n=length(tag_ini); % Number of enclosed tag structures
        % Extract tag names and contents
        if n>0
            for i=1:n % Loop through all tags
                tag_close = (XML(tag_ini(i)+1)=='/'); % 1 for closing contents and 0 for opening
                tmp_tag = XML(tag_ini(i)+1+tag_close:tag_end(i)-1);
                tag_contents{i,1} = tmp_tag;  % First column contains names
            end
        end
        
        % Extract strings containing start time information
        t_start_idx = strncmp('START STATE',tag_contents,11);
        t_start_idx = find(t_start_idx == 1)-1; % Get start event time index
        str_start = tag_contents{t_start_idx,1}; % String containing start datevec info
        
        % Extract datevec start from the start event string
        tmpIndexes = [strfind(str_start,'"'),strfind(str_start,',')];
        tmpIndexes = sort(tmpIndexes);
        datevec_start = zeros(1,6);
        for i=1:length(tmpIndexes)-1
            tmpI = (tmpIndexes(i)+1):(tmpIndexes(i+1)-1);
            datevec_start(1,i) = str2double(str_start(tmpI));
        end
        
        % Return .xml file in MATLAB datenum format
        output = datenum(datevec_start);
    end

    function s=file2str(x)
        %FILE2STR reads textfile into a single long string
        %
        % Syntax: s=file2str(x)
        %
        % Description
        %   x is a filename
        %   s is the long string with all contents
        %
        % Jonas Almeida, almeidaj@musc.edu, 30 June 2003, MAT4NAT Tbox
        % Fixed to return an empty string if the file is empty
        % mj July 2012
        
        s = '' ;
        fid=fopen(x,'r');
        if fid<0, return, end
        i=1;
        y = {} ;
        while ~feof(fid)
            yy = fgetl(fid);
            if ischar(yy)
                y{i} = yy ;
                i=i+1;
            end
        end
        fclose(fid);
        if length(y)>0
            s=strcat(y{:});
        end
    end

end