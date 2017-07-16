function dtagsignal()
% MS 2016.07.11
% Custom graphic user interface designed for reviewing calls marked with
% the dtagmark.m interface, especially optimized for low-frequency blue
% whale calls. The dtagsignal.m GUI allows for signals to be precisely
% selected and saved in a .mat file for easy parameter calculations.

clear; close all;

% Initialize data structures for deployment data
AUD = []; % Initialize structure for audio deployment data
ACC = []; % Initialize structure for accel deployment data
ACC_raw = ACC; % Initialize structure for unfiltered accel data

% Initialize data structures for the selected signal
AUDIOSIGNAL = struct('y',[],'T',[],'fs',[],...
    'noise_pre',[],'noise_post',[],...
    'y_raw',[],'T_raw',[],...
    'file',[],'call',[],'comment',[]);
ACCELSIGNAL = struct('y',[],'T',[],'fs',[],...
    'noise_pre',[],'noise_post',[],...
    'y_raw',[],'T_raw',[],...
    'bandpass',[],'p',[]);

% Initialize data structures for current window of deployment data
WINDOW = struct('time',[],'dur',[],'rect',[],...
    'file',[],'call',[],'comment',[]);

% Initialize output cell array and define signal selection parameters
savedSignals = []; % Output cell array for saved signals {AUDIO,ACCEL}
counter = 1; % Count how many signals have been selected
signalList = []; % Short form of savedSignals to display on GUI
outputFile = 'G:\Mark\data_bm\TEMP.mat'; % Out file to save signal array   <--- PROVIDE OUTPUT FILENAME FOR SIGNALS!!
noiseDuration = 1; % Duration of noise sample pre and post signal [s]    <--- PROVIDE DURATION OF NOISE SAMPLES!!

% Fourier analysis parameters
aud_pt = 512; % FFT size for audio spectrogram
aud_ovl = 0.90; % Fraction of overlap for audio spectrogram
acc_pt = 128; % FFT size for accel spectrogram
acc_ovl = 0.99; % Fraction of overlap for audio spectrogram

% Accelerometer data filter parameters
bandpass = [0,0]; % Accelerometer high-pass filter cutoff frequency (Hz)
filtOrder = 6;

% Audit log variables and display parameters
log_filename = []; % Audit log file name
log = {'DATENUM','DURATION','FILE','DATESTR','CALL','COMMENT'};
auditList = []; % Short form of log entries to display on GUI
commLength = 16; % Number of comment characters to display in log entries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Graphic User Interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up figure and panels
fig = figure('Name','DTAGsignal','Position',[50,100,1200,800]);
pan_c = uipanel('Parent',fig,'Position',[.01 .01 .18 .98],'Title','Controls');
pan_p = uipanel('Parent',fig,'Position',[.20 .01 .79 .98],'Title','Plot');
pan_c_filt = uipanel('Parent',pan_c,'Position',[.01 .42 .98 .12],...
    'Title','Accelerometer band-pass filter');

% Build GUI control buttons
c_loaddeployment = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .96 .96 .03],...
    'String','Load Deployment Data','Callback',@callback_pushbutton);
c_loadaudit = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .92 .96 .03],...
    'String','Load Audit Log','Callback',@callback_pushbutton);
c_auditlist = uicontrol('Parent',pan_c,'Style','listbox',...
    'Units','Normalized','Position',[.02 .56 .96 .35],...
    'String',auditList,'Callback',@callback_listbox);
c_addsignal = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .34 .96 .06],...
    'String','Add signal to list','Callback',@callback_pushbutton);
c_signallist = uicontrol('Parent',pan_c,'Style','listbox',...
    'Units','Normalized','Position',[.02 .08 .96 .25],...
    'String',signalList,'Callback',@callback_listbox);
c_savelist = uicontrol('Parent',pan_c,'Style','pushbutton',...
    'Units','Normalized','Position',[.02 .01 .96 .06],...
    'String',['Save signals: ',outputFile],'Callback',@callback_pushbutton);
% Accelerometer band-pass filter parameters
uicontrol('Parent',pan_c_filt,'Style','Text',...
    'Units','Normalized','Position',[.01 .52 .59 .30],...
    'String','Frequency band cutoffs (Hz)');
c_bandpass1 = uicontrol('Parent',pan_c_filt,'Style','edit',...
    'Units','Normalized','Position',[.60 .56 .19 .30],...
    'String',num2str(bandpass(1)),'Callback',@callback_edit);
c_bandpass2 = uicontrol('Parent',pan_c_filt,'Style','edit',...
    'Units','Normalized','Position',[.79 .56 .19 .30],...
    'String',num2str(bandpass(2)),'Callback',@callback_edit);
c_filtrun = uicontrol('Parent',pan_c_filt,'Style','pushbutton',...
    'Units','Normalized','Position',[.01 .08 .76 .30],...
    'String','Butterworth filter - order:','Callback',@callback_pushbutton);
c_filtorder = uicontrol('Parent',pan_c_filt,'Style','edit',...
    'Units','Normalized','Position',[.80 .08 .19 .30],...
    'String',num2str(filtOrder),'Callback',@callback_edit);

% Set up plot panel format
p_header = uicontrol('Parent',pan_p,'Style','Text',...
    'Units','Normalized','Position',[.01 .96 .98 .03],...
    'String',log_filename);
% Set up plot panel axes
ax_audio = gobjects(1,2);
ax_accel = gobjects(3,2);
ax_audio(1,1) = axes('Parent',pan_p,'Position',[.05 .57 .43 .38]);
ax_audio(1,2) = axes('Parent',pan_p,'Position',[.55 .57 .43 .38]);
ax_accel(1,1) = axes('Parent',pan_p,'Position',[.05 .05 .43 .15]);
ax_accel(2,1) = axes('Parent',pan_p,'Position',[.05 .20 .43 .15]);
ax_accel(3,1) = axes('Parent',pan_p,'Position',[.05 .35 .43 .15]);
ax_accel(1,2) = axes('Parent',pan_p,'Position',[.55 .05 .43 .15]);
ax_accel(2,2) = axes('Parent',pan_p,'Position',[.55 .20 .43 .15]);
ax_accel(3,2) = axes('Parent',pan_p,'Position',[.55 .35 .43 .15]);
% Link axes to align signals across plots
linkaxes([ax_audio;ax_accel],'x')
linkaxes(ax_accel(:,1),'y')
linkaxes(ax_accel(:,2),'y')
% Set buttondown function of audio plots to select signal for calculations
set(ax_audio,'ButtonDownFcn',@buttondown_audio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI Callback Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function callback_pushbutton(source,~)
        switch source
            case c_loaddeployment
                % Request user to select directory with .wav and .xml files
                recDir = uigetdir('G:\Mark\data_bm\',...
                    'Select directory with .wav and .xml files');
                if recDir == 0
                    disp('Load operation canceled'); return
                end
                recDir = [recDir,'\'];
                
                % Call DATA HANDLING function load_deployment to read-in
                % audio and accelerometer data
                [AUD, ACC] = load_deployment(recDir);
                ACC_raw = ACC; % Store raw (unfiltered) accelerometer data
            
            case c_loadaudit
                % Request user to open log audit (.xlsx) file
                [FileName,PathName] = uigetfile('G:\Mark\data_bm\*.xlsx',...
                    'Select audit log (Excel spreadsheet)');
                log_filename = [PathName,FileName];
                if ~ischar(log_filename)
                    disp('Load operation canceled'); return
                end
                
                disp(['Loading audit log file: ',log_filename])
                [num,txt,raw] = xlsread(log_filename);
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
                
                % Display audit log file name on plot panel header
                set(p_header,'String',log_filename)
                
            case c_filtrun
                % If cutoff frequencies are provided and accelerometer data
                % has been loaded, apply filter to raw accelerometer data
                if isempty(ACC)
                   warning('No accelerometer data found'); return 
                end
                
                if 0<bandpass(1) && bandpass(1)<ACC.fs
                    if bandpass(2) <= bandpass(1)
                        % High-pass filter if only low cutoff provided
                        [b,a] = butter(filtOrder,bandpass(1)/(ACC.fs/2),'high');
                        for i=1:length(ACC.y)
                            ACC.y{i} = filter(b,a,ACC_raw.y{i});
                        end
                        disp('Acceleromter data high-pass filtered')
                    elseif bandpass(1)<bandpass(2) && bandpass(2)<ACC.fs
                        % Band-pass filter if both cutoffs provided
                        [b,a] = butter(filtOrder,bandpass/(ACC.fs/2),'bandpass');
                         for i=1:length(ACC.y)
                            ACC.y{i} = filter(b,a,ACC_raw.y{i});
                         end
                         disp('Acceleromter data band-pass filtered')
                    end
                else % If no valid cutoffs, restore unfilterd data
                    ACC.y = ACC_raw.y;
                    disp('Acceleromter data not filtered')
                end
                
                % If a window is already displayed, refresh plots
                if isempty(WINDOW.time)
                    return % return if no window currently diplayed
                end
                time = WINDOW.time;
                dur = WINDOW.dur;
                % Call DATA HANDLING function locate_signal to find marked
                % signal in the audio and acceleromter data
                [audioChip,audioSignal] = locate_signal(AUD.T,time,dur);
                [accelChip,accelSignal] = locate_signal(ACC.T,time,dur);
                % Call DATA PLOTTING function plot_signal to visualize
                % signal waveforms and spectrograms
                if ~isempty(audioSignal)
                    plot_signal(ax_audio,AUD.T{audioChip}(audioSignal),...
                        AUD.y{audioChip}(audioSignal,:),...
                        AUD.fs,aud_pt,aud_ovl)
                    set(ax_audio(1,1),'ButtonDownFcn',@buttondown_audio);
                    set(imhandles(ax_audio(1,2)),'ButtonDownFcn',@buttondown_audio);
                end
                if ~isempty(accelSignal)
                    plot_signal(ax_accel,ACC.T{accelChip}(accelSignal),...
                        ACC.y{accelChip}(accelSignal,:),...
                        ACC.fs,acc_pt,acc_ovl)
                end
                
            case c_addsignal
                % If selection is valid, add signal to the savedSignals
                % cell array and update the GUI signalList
                if length(AUDIOSIGNAL.T) > 1
                    savedSignals{counter,1} = AUDIOSIGNAL;
                    savedSignals{counter,2} = ACCELSIGNAL;
                    counter = counter + 1;
                    
                    % Add signal to the GUI signalList listbox
                    timeString = datestr(AUDIOSIGNAL.T(1),...
                        'yyyy/mm/dd HH:MM:SS');
                    callString = [AUDIOSIGNAL.call,'      '];
                    callString = callString(1:6);
                    signalList = [signalList;timeString,'   ',callString];
                    set(c_signallist,'String',signalList)
                    
                    % Call the callback function for the newly added line
                    % to the listbox (to allow user to edit call / comment)
                    set(c_signallist,'Value',size(signalList,1))
                    callback_listbox(c_signallist)
                    
                    % Use color to prompt user to save signal list
                    set(c_addsignal,'BackgroundColor',[.94 .94 .94])
                    set(c_savelist,'BackgroundColor',[.94 .64 .64])
                else
                    warning('No valid signal selected')
                end
                
            case c_savelist
                % Save the output cell array savedSignals to the output
                % .mat file specified at the top of this script
                save(outputFile,'savedSignals','AUD','ACC')
                set(c_savelist,'BackgroundColor',[.94 .94 .94])
        end
    end

 function callback_edit(source,~)
        switch source
            case c_bandpass1 % Low frequency cutoff
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    bandpass(1) = abs(val);
                end
            case c_bandpass2 % High frequency cutoff
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    bandpass(2) = abs(val);
                end
            case c_filtorder
                val = str2double(get(source,'String'));
                if ~isnan(val)
                    filtOrder = floor(abs(val));
                end
        end
    end

    function callback_listbox(source,~)
        switch source
            case c_auditlist
                % Select window to view from the audit log
                if isempty(AUD)
                   warning('No deployment data loaded'); return
                end
                entryIdx = get(source,'Value');
                if isempty(entryIdx)
                    warning('No audit log loaded'); return
                else
                    time = log{entryIdx+1,1};
                    dur = log{entryIdx+1,2};
                end
                
                % Store current audit log entry information in WINDOW
                WINDOW.time = time;
                WINDOW.dur = dur;
                WINDOW.file = log{entryIdx+1,3};
                WINDOW.call = log{entryIdx+1,5};
                WINDOW.comment = log{entryIdx+1,6};
                
                % Call DATA HANDLING function locate_signal to find marked
                % signal in the audio and acceleromter data
                [audioChip,audioSignal] = locate_signal(AUD.T,time,dur);
                [accelChip,accelSignal] = locate_signal(ACC.T,time,dur);
                
                % Call DATA PLOTTING function plot_signal to visualize
                % signal waveforms and spectrograms
                if ~isempty(audioSignal)
                    plot_signal(ax_audio,AUD.T{audioChip}(audioSignal),...
                        AUD.y{audioChip}(audioSignal,:),...
                        AUD.fs,aud_pt,aud_ovl)
                    set(ax_audio(1,1),'ButtonDownFcn',@buttondown_audio);
                    set(imhandles(ax_audio(1,2)),'ButtonDownFcn',@buttondown_audio);
                end
                if ~isempty(accelSignal)
                    plot_signal(ax_accel,ACC.T{accelChip}(accelSignal),...
                        ACC.y{accelChip}(accelSignal,:),...
                        ACC.fs,acc_pt,acc_ovl)
                end
                
            case c_signallist
                % Edit the call type and comment of a saved signal
                entryIdx = get(source,'Value');
                if isempty(entryIdx)
                    warning('No signal selected'); return
                else
                    % If the signal is in the current window, re-plot the
                    % rectangle lines used to select signal
                    if WINDOW.time <= savedSignals{entryIdx,1}.T(1) && ...
                            savedSignals{entryIdx,1}.T(1) <= ...
                            WINDOW.time + WINDOW.dur
                        
                        tmpT = (savedSignals{entryIdx,1}.T(1)...
                            - WINDOW.time)*24*60*60;
                        tmpDur = (savedSignals{entryIdx,1}.T(end)...
                            - savedSignals{entryIdx,1}.T(1))*24*60*60;
                        tmpRect = [tmpT,1,tmpDur,1];
                        
                        plot_rectlines(ax_audio(1,1),tmpRect)
                        plot_rectlines(ax_audio(1,2),tmpRect)
                        plot_rectlines(ax_accel(1,1),tmpRect)
                        plot_rectlines(ax_accel(2,1),tmpRect)
                        plot_rectlines(ax_accel(3,1),tmpRect)
                        plot_rectlines(ax_accel(1,2),tmpRect)
                        plot_rectlines(ax_accel(2,2),tmpRect)
                        plot_rectlines(ax_accel(3,2),tmpRect)
                    end
                    
                    % Prompt user to edit the call type and comment
                    prompt = {'Call Type:','Comment:'};
                    predef = {savedSignals{entryIdx,1}.call,...
                        savedSignals{entryIdx,1}.comment};
                    response = inputdlg(prompt,'Edit signal',1,predef);
                    if ~isempty(response)
                        savedSignals{entryIdx,1}.call = response{1};
                        savedSignals{entryIdx,1}.comment = response{2};
                    end
                    timeString = datestr(savedSignals{entryIdx,1}.T(1),...
                        'yyyy/mm/dd HH:MM:SS');
                    callString = [savedSignals{entryIdx,1}.call,'      '];
                    callString = callString(1:6);
                    signalList(entryIdx,:) = [timeString,'   ',callString];
                    set(c_signallist,'String',signalList)
                end
        end
    end

    function buttondown_audio(source,event)
        % Upon right-clicking the axes, user is able to draw a rectangle to
        % select the call of interest within the window. The time vector
        % and audio signal along with the rectangle dimensions will be
        % stored in the data structure 'SELECT'. The signal stored in
        % 'SELECT' may be used for further calculations
        if event.Button == 3
            % If source is not an axes object, get the axes object
            if ~strcmp(get(source,'Type'),'axes')
                source = get(source,'Parent');
            end
            % Request new rectangle
            delete(findobj(get(source,'Children'),'Type','Rectangle'));
            rect = getrect(source);
            % Store selection data and plot rectangle horizontal lines
            if ~isempty(WINDOW.time)
                WINDOW.rect = rect;
                time = WINDOW.time + rect(1)/(24*60*60);
                dur = rect(3)/(24*60*60);
                
                % Store selected signal data in SIGNAL structures
                [audioChip,audioSignal] = locate_signal(AUD.T,time,dur);
                AUDIOSIGNAL.y = AUD.y{audioChip}(audioSignal,:);
                AUDIOSIGNAL.T = AUD.T{audioChip}(audioSignal);
                AUDIOSIGNAL.fs = AUD.fs;
                AUDIOSIGNAL.file = WINDOW.file;
                AUDIOSIGNAL.call = WINDOW.call;
                AUDIOSIGNAL.comment = WINDOW.comment;
                if ~isempty(ACC.fs)
                    [accelChip,accelSignal] = locate_signal(ACC.T,time,dur);
                    ACCELSIGNAL.y = ACC.y{accelChip}(accelSignal,:);
                    ACCELSIGNAL.T = ACC.T{accelChip}(accelSignal);
                    ACCELSIGNAL.fs = ACC.fs;
                    ACCELSIGNAL.bandpass = bandpass;
                    ACCELSIGNAL.p = ACC.p{accelChip}(accelSignal,:);
                end
                
                % Collect noise sample before and after selected signal
                noiseDur = noiseDuration/(24*60*60);
                preTime = AUDIOSIGNAL.T(1) - noiseDur;
                postTime = AUDIOSIGNAL.T(end);
                [preC,preS] = locate_signal(AUD.T,preTime,noiseDur);
                [postC,postS] = locate_signal(AUD.T,postTime,noiseDur);
                AUDIOSIGNAL.noise_pre = AUD.y{preC}(preS,:);
                AUDIOSIGNAL.noise_post = AUD.y{postC}(postS,:);
                if ~isempty(ACC.fs)
                    preTime = ACCELSIGNAL.T(1) - noiseDur;
                    postTime = ACCELSIGNAL.T(end);
                    [preC,preS] = locate_signal(ACC.T,preTime,noiseDur);
                    [postC,postS] = locate_signal(ACC.T,postTime,noiseDur);
                    ACCELSIGNAL.noise_pre = ACC.y{preC}(preS,:);
                    ACCELSIGNAL.noise_post = ACC.y{postC}(postS,:);
                end
                
                % Collect raw signal from entire window
                [windowAudioChip,windowAudioSignal] = locate_signal(AUD.T,WINDOW.time,WINDOW.dur);
                [windowAccelChip,windowAccelSignal] = locate_signal(ACC.T,WINDOW.time,WINDOW.dur);
                AUDIOSIGNAL.y_raw = AUD.y{windowAudioChip}(windowAudioSignal);
                AUDIOSIGNAL.T_raw = AUD.T{windowAudioChip}(windowAudioSignal);
                ACCELSIGNAL.y_raw = ACC_raw.y{windowAccelChip}(windowAccelSignal,:);
                ACCELSIGNAL.T_raw = ACC.T{windowAccelChip}(windowAccelSignal);
                
                % Plot rectangle horizontal lines on desired plots
                plot_rectlines(ax_audio(1,1),WINDOW.rect)
                plot_rectlines(ax_audio(1,2),WINDOW.rect)
                plot_rectlines(ax_accel(1,1),WINDOW.rect)
                plot_rectlines(ax_accel(2,1),WINDOW.rect)
                plot_rectlines(ax_accel(3,1),WINDOW.rect)
                plot_rectlines(ax_accel(1,2),WINDOW.rect)
                plot_rectlines(ax_accel(2,2),WINDOW.rect)
                plot_rectlines(ax_accel(3,2),WINDOW.rect)
                
                % Use color to prompt user to save the selected signal
                set(c_addsignal,'BackgroundColor',[.64 .94 .64])
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_signal(ax,t,y,fs,pt,ovl)
        % Function plots signal waveforms and spectrograms next to each
        % other on the provided axes.
        
        % Convert time vector into seconds
        t = (t-t(1))*24*60*60;
        
        for row = 1:size(ax,1)
            % Plot waveforms in the left column
            plot(ax(row,1),t,y(:,row))
            if row ~= 1
                set(ax(row,1),'xtick',[]) % Hide x-axis labels
            end
            % Calculate spectrograms in the right column
            [S,F,T] = spectrogram(y(:,row),...
                hamming(pt),floor(ovl*pt),pt,fs);
            S = 10*log10(abs(S)); % Convert to decibels
            % Plot spectrogram and format axis
            imagesc(T,F,S,'Parent',ax(row,2))
            colormap(ax(row,2),'bone')
            set(ax(row,2),'YDir','Normal')
            if row ~= 1
                set(ax(row,2),'xtick',[]) % Hide x-axis labels
            end
        end
    end

    function plot_rectlines(ax,rect)
        % Plot full-height horiztonal lines of the provided rectangle on
        % the specified axes
        delete(findobj(get(ax,'Children'),'Type','Rectangle'));
        
        ylim = get(ax,'YLim');
        rect(2) = ylim(1);
        rect(4) = ylim(2)- ylim(1);
        
        hold(ax,'on')
        rectangle('Parent',ax,'Position',rect,'EdgeColor',[1 0 0]);
        hold(ax,'off')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Handling Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [AUDIO,ACCEL] = load_deployment(recDir)
        % Function will return structures containing the aligned audio and
        % accelerometer data for the deployment indicated by the input
        % directory. The input directory (recDir) should point to the
        % directory containing a deployment's .wav and .xml files.
        
        % Output structures (AUDIO and ACCEL) each contain three fields:
        %   AUDIO.y = signal vectors (cell array)
        %   AUDIO.T = time vectors (cell array) [MATLAB datenums]
        %   AUDIO.fs = sample rate [Hz]
        AUDIO = struct('y',[],'T',[],'fs',[]);
        ACCEL = struct('y',[],'T',[],'fs',[],'p',[]);
        
        % Get names of .wav and .xml files in deployment
        wavFiles = dir(fullfile(recDir,'*.wav'));
        xmlFiles = dir(fullfile(recDir,'*.xml'));
        if isempty(wavFiles) || isempty(xmlFiles)
            warning('.wav and .xml files not found')
            return
        end
        
        % Get start times of all .xml times in deployment
        xmlTimes = zeros(length(xmlFiles),1);
        for itr0 = 1:length(xmlTimes)
            xmlFileTmp = [recDir,xmlFiles(itr0).name]; % .xml file name
            xmlTimes(itr0,1) = extractxmlstart(xmlFileTmp); % local fcn
        end
        
        % Initialize cell arrays for audio signal and time vectors
        y = cell(length(wavFiles),1);
        T = cell(length(wavFiles),1);
        fs = 0;
        
        % Read-in audio signals from .wav files
        for itr0 = 1:length(wavFiles)
            tmpFile = [recDir,wavFiles(itr0).name];
            disp(['Reading: ',tmpFile, ' (file ',num2str(itr0),...
                ' of ',num2str(length(wavFiles)),')'])
            [y{itr0},fs] = audioread(tmpFile);
            % Set up time vectors for audio files (MATLAB datenums)
            T{itr0} = 1:(length(y{itr0}));
            T{itr0} = (T{itr0}/fs)/(60*60*24) + xmlTimes(itr0);
        end
        
        % Store audio data in AUDIO structure
        AUDIO.fs = fs;
        AUDIO.T = T;
        AUDIO.y = y;
        
        % Find deployment's prh.mat file with accelerometer data
        slashIdx = strfind(recDir,'\');
        prhFile = [recDir(1:slashIdx(end-2)),'prh',...
            recDir(slashIdx(end-1):end-1),'prh.mat'];
        if exist(prhFile,'file') == 0
            warning('prh.mat file not found')
            return
        end
        
        % Initialize cell arrays for accelerometer signals and time vectors
        Ay = cell(length(wavFiles),1);
        AT = cell(length(wavFiles),1);
        Ap = cell(length(wavFiles),1);
        
        % Read-in accelerometer data from prh.mat file
        disp(['Loading PRH file data: ',prhFile])
        m = matfile(prhFile); % Load acc data, fs, and depth
        Ay_tmp = m.A; Afs = m.fs; Ap_tmp = m.p;
        % (1) Set up time vectors for accelerometer data (T_a) to
        % align with the time vectors for audio data (T)
        % (2) Set up cell array for accelerometer data (A) to align
        % with the audio data (y)
        for itr0 = 1:length(wavFiles)
            % Align accelerometer time vectors with audio time vectors
            t_start = xmlTimes(itr0);
            t_end = T{itr0}(end);
            AT{itr0} = t_start:(1/Afs)/(60*60*24):t_end;
            % Align accelerometer signal vectors with audio signal vectors
            A_end = min(length(Ay_tmp),length(AT{itr0}));
            Ay{itr0,1} = Ay_tmp(1:A_end,:);
            Ay_tmp = Ay_tmp((length(AT{itr0})+1):end,:);
            
            Ap{itr0,1} = Ap_tmp(1:A_end,:);
            Ap_tmp = Ap_tmp((length(AT{itr0})+1):end,:);
        end
        
        % Store accelerometer data in ACCEL structure
        ACCEL.fs = Afs;
        ACCEL.T = AT;
        ACCEL.y = Ay;
        ACCEL.p = Ap;
    end

    function [chipIdx,signalIdx] = locate_signal(T,time,dur)
        % Function will return the chip index and the vector of time
        % indexes (within input 'T') corresponding to the signal specified
        % by 'time' and 'dur'. Inputs:
        %   T = cell array of time vectors [MATLAB datenums]
        %   time = signal start time [MATLAB datenum]
        %   dur = signal duration [MATLAB datenum]
        if isempty(T)
           chipIdx = []; signalIdx = []; return 
        end
        
        chipStartTimes = zeros(1,length(T));
        for i = 1:length(T)
            chipStartTimes(i) = T{i}(1);
        end
        chipIdx = find(chipStartTimes <= time,1,'last');
        signalStart = find(T{chipIdx} <= time,1,'last');
        signalEnd = find(T{chipIdx} <= time+dur,1,'last');
        signalIdx = signalStart:signalEnd;
    end

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

    function s = file2str(x)
        % FILE2STR reads textfile into a single long string
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
            if ischar(yy),
                y{i} = yy ;
                i=i+1;
            end
        end
        fclose(fid);
        if length(y)>0,
            s=strcat(y{:});
        end
    end

end