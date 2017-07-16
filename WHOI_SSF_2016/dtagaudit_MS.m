function [RES,AXs,hhh] = dtagaudit_MS(tag,tcue,RES)
% MS 2016.06.14
%
%     R = dtagaudit_synch(tag,tcue,R)
%
%     DTAG audit tool compatible with DTAG versions 1 through 3
%     [Support for simultaneous tags plus multiple smaller hacks]
%
%     tag is the tag deployment string e.g., 'sw03_207a'
%     tcue is the time in seconds-since-tag-on to start displaying from
%     R is an optional audit structure to edit or augment
%     Output:
%        R is the audit structure made in the session. Use saveaudit
%        to save this to a file.
%
%     OPERATION
%     Type or click on the display for the following functions:
%     - type 'f' to go to the next block
%     - type 'b' to go to the previous block
%     - click on the graph to get the time cue, depth, time-to-lastf
%       and frequency of an event. Time-to-last is the elapsed time
%       between the current click point and the point last clicked.
%       Results display in the matlab command window.
%     - type 's' to select the current segment and add it to the audit.
%       You will be prompted to enter a sound type on the matlab command
%       window. Enter a single word and type return when complete.
%     - type 'l' to select the currect cursor position and add it to the
%       audit as a 0-length event. You will be prompted to enter a sound
%       type on the matlab command window. Enter a single word and type
%       return when complete.
%     - type 'x' to change the audit entry at the cursor position.
%       If there is no audit entry at the cursor, nothing happens.
%       Otherwise, a prompt will appear with the current cue type
%       for all cues around cursor position. Delete text to remove that
%       cue, or change it to update it
%     - type 'i' to play the current segment and simultaneously save to
%       current folder as tempsound.wav. Rename to keep sound file.
%     - type 'p' to play the displayed sound segment
%       through the computer speaker/headphone jack.
%     - type 'q' or press the right hand mouse button to finish auditing.
%     - type 'a' to report the angle of arrival of the selected segment
%
%       Multiple tag comparison package
%     - type 'o' to adjust tagaudit_synch parameters (tag ID and time shift)
%     - type 't' to make a synched tagaudit plot of the current segment
%     - type 'T' to make a synched tagaudit angle-of-arrival plot
%     - type 'S' to mark non-focal sounds
%
%
%     mark johnson, WHOI
%     majohnson@whoi.edu
%     last modified March 2005
%     added buttons and updated audit structure
%
%     Added by Frants Jensen, AU and WHOI:
%     - type 'i' to play sound only in the marked interval,
%       and simultaneously save the interval to a temporary
%       wave-file named 'tempsound.wav' (rename to use in presentations)
%     - changed function 'x' to change/remove any cue type near cursor
%
%     Added option for choosing maximum frequency to display in specgram
%
%     Added compatibility with 4-channel DTAG files
%
%     Multiple tag comparison package (separate files required):
%     - type 'o' to set up multiple tag comparison.
%     - type 't' to compare spectrograms and received levels across tags
%       (will open separate window)
%     - type 'T' to compare angle-of-arrival across tags (separate window)
%     - type 'S' to start marking sequences on non-focal tags (requires 't' first):
%                Press twice in same axis to activate sequence
%                Press 's' to save active sequence as cue
%                Press 'l' to mark instant point
%                Press 'x' to edit any cue types close to cursor
%                   -change text to update cue type
%                   -remove text to delete cue
%                Press any other button to go back to auditing window
%                WARNING: This will load, change, and save audits!
%                Make sure to have newest audits when making changes
%     Multiple tag comparison package requires:
%     d3tagaudit_synch_int.m
%     d3tagaudit_synch_aoa_int.m
%     d3resamplesound.m

NS = 60;          % number of seconds to display (60)
BL = 1*512 ;      % specgram (fft) block size (4*512)
CLIM = [-70 -10] ; % color axis limits in dB for specgram ([-90 -10])
CH = 1 ;           % which channel to display if multichannel audio
THRESH = 0 ;       % click detector threshold, 0 to disable
volume = 8 ;       % amplification factor for audio output - often needed to
% hear weak signals (if volume>1, loud transients will
% be clipped when playing the sound cut
SOUND_FH = 0 ;   % high-pass filter for sound playback - 0 for no filter (200)
SOUND_DF = 10 ;     % decimation factor for playing sound
MAXYONCLICKDISPLAY = 0.1 ;
MAXYONFREQDISPLAY  = 1 ;    % In kHz (but max half of AFS_RESAMPLE)
AFS_RESAMPLE       = 2000 ; % Resample sound to limit data and speed up specgram (2000)
QUERY_LABEL        = 0;      % Ask for extra label when pressing "i" and saving sound

% Default Settings for multiple tag comparisons
possibletags = '' ; othertags    = '' ; AXS_synch    = [] ; time_shift   = 0  ;

% Find tag version
tagver = dtagtype(tag);

% Decide on hydrophone separation
if tagver==3,
    AOA_SCF = 1500 / (0.045) ; % sound speed / hydrophone separation (adjusted for D3)
elseif tagver==2,
    AOA_SCF = 1500/0.0228 ; % v/h % This should be /0.0228 (0.9 inch)
end

% Check if there's a recent synch audit saved settings:
synch_fname = [tag(1:9) '_synch.mat'] ;
if exist(synch_fname),
    load(synch_fname)
    disp(['Loading settings from ' synch_fname])
end

% high-pass filter frequencies (Hz) for click detector
switch tag(1:2),
    case 'zc',      % for ziphius use:
        FH = 20000 ;            % high-pass filter settings for click detection
        TC = 0.5e-3 ;           % power averaging time constant in seconds
    case 'md',      % for mesoplodon use:
        FH = 20000 ;
        TC = 0.5e-3 ;
    case 'pw',      % for pilot whale use:
        FH = 10000 ;
        TC = 0.5e-3 ;
    case 'gm',      % for long-finned pilot whale use:
        FH = 10000 ;
        TC = 0.5e-3 ;
    case 'sw',      % for sperm whale use:
        FH = 3000 ;
        TC = 2.5e-3 ;
    case 'bm'       % for blue whale use: <--- MS's edit
        FH = 5 ;
        TC = 2.5e-3 ;
    otherwise,      % for others use:
        FH = 5000 ;
        TC = 0.5e-3 ;
end

% Check if audit structure exists
if nargin<3 | isempty(RES),
    RES.cue = [] ;
    RES.comment = [] ;
end

%%%%%%%%%%%%%%%%%%% Load prh file %%%%%%%%%%%%%%%%%%%%%
k = loadprh(tag,'p','fs') ;           % read p and fs from the sensor file
if k==0,
    fprintf('Unable to find a PRH file - continuing without\n') ;
    p = [] ; fs = [] ;
end


%%%%%%%%%%%%%%%%%%% Prepare filters %%%%%%%%%%%%%%%%%%%%%

% check sampling rate and reset to AFS_RESAMPLE if too high
[x,afs] = dtagwavread(tag,tcue,0.1); % dtagwavread(tag,tcue,0.1);
if afs>AFS_RESAMPLE, afs=AFS_RESAMPLE; end

% filters for sound playback
if SOUND_FH > 0,
    [bs as] = butter(6,SOUND_FH/(afs/2),'high') ;
else
    bs = [] ;
end

% high pass filter and smoothing filter for envelope
% Chebyshev Type I filter:
% cheby1(Filter order, Peak-to-peak passband ripple, Passband edge frequency, Filter type)
[bh ah] = cheby1(6,0.5,FH/(afs/2),'high') ;
pp = 1/TC/afs ;


%%%%%%%%%%%%%%%%% Prepare figure plot %%%%%%%%%%%%%%%%%%%

current = [0 0] ;
figure(1),clf
if ~isempty(p),
    kb = 1:floor(NS*fs) ;
    AXm = axes('position',[0.11,0.76,0.78,0.18]) ;
    AXc = axes('position',[0.11,0.70,0.78,0.05]) ;
    AXs = axes('position',[0.11,0.34,0.78,0.35]) ;
    AXp = axes('position',[0.11,0.11,0.78,0.2]) ;
else
    AXm = axes('position',[0.11,0.60,0.78,0.34]) ;
    AXc = axes('position',[0.11,0.52,0.78,0.07]) ;
    AXs = axes('position',[0.11,0.11,0.78,0.38]) ;
end

bc = get(gcf,'Color') ;
set(AXc,'XLim',[0 1],'YLim',[0 1]) ;
set(AXc,'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;

while 1,
    [x,afs_org] = dtagwavread(tag,tcue,NS);
    if isempty(x), return, end
    
    % Else, resample to conserve memory
    [x,afs] = d3resamplesound(x,afs_org,CH,AFS_RESAMPLE) ;
    
    % Construct spectrogram - use top line for Matlab version w/o specgram.m
    % [B F T] = spectrogram(x,hamming(BL),floor(BL/1.3),BL,afs, 'yaxis') ;
    [B F T] = specgram(x,BL,afs,hamming(BL),floor(BL/1.001)) ;
    % [B F T] = spectrogram(x,hamming(BL),floor(BL/1.3)) ;
    
    % If a click detection threshold is defined, run click detector
    if THRESH,
        [cl xx] = findclicks(x,THRESH,afs,FH) ;
    end
    
    % Create top plot
    xx = filter(pp,[1 -(1-pp)],abs(filter(bh,ah,x))) ;
    kk = 1:5:length(xx) ;
    axes(AXm), plot(tcue+kk/afs,10*log10(xx(kk)),'k') ; grid
    set(AXm,'XAxisLocation','top') ;
    axis([tcue tcue+NS -40 10*log10(MAXYONCLICKDISPLAY)]) ;
    ylabel('Intensity, dB')
    
    % Add audit cues
    plotRES(AXc,RES,[tcue tcue+NS]) ;
    
    % Create bottom plot
    if ~isempty(p),
        ks = kb + round(tcue*fs) ;
        axes(AXp),plot(ks/fs,p(ks)), grid
        set(gca,'YDir','reverse') ;
        axis([tcue tcue+max(T) get(gca,'YLim')]) ;
        xlabel('Time, s'), ylabel('Depth, m')
    end
    
    % Construct spectrogram plot
    BB = adjust2Axis(20*log10(abs(B))) ;
    axes(AXs), imagesc(tcue+T,F/1000,BB,CLIM) ; axis xy, grid ;
    colormap('winter')
    if ~isempty(p),
        set(AXs,'XTickLabel',[]) ;
    else
        xlabel('Time, s')
    end
    ylabel('Frequency, kHz')
    
    % Limit frequency axis on spectrogram
    yl = get(gca,'YLim') ;
    yl(2) = min([yl(2) MAXYONFREQDISPLAY]) ;
    axis([tcue tcue+NS yl]) ;
    
    % Add current selection to spectrogram plot
    hold on
    hhh = plot([0 0],0.8*min([afs/2000 MAXYONFREQDISPLAY])*[1 1],'k*-') ;    % plot cursor
    if THRESH,
        plot(cl+tcue,0.4*afs*ones(length(cl),1),'wo') ;
    end
    hold off
    
    
    done = 0 ;
    while done == 0,
        axes(AXs) ; pause(0) ;
        [gx gy button] = ginput(1) ;
        
        % QUIT AUDITING PROGRAM AND RETURN CUE STRUCTURE
        if button==3 | button=='q',
            save tagaudit_RECOVER RES
            disp('Updating R structure after tagaudit. Remember saveaudit(tag,R)')
            return
            
            % INSERT SEQUENCE CUE
        elseif button=='s',
            ss = input(' Enter comment... ','s') ;
            cc = sort(current) ;
            RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
            RES.stype{size(RES.cue,1)} = ss ;
            save tagaudit_RECOVER RES
            plotRES(AXc,RES,[tcue tcue+NS]) ;
            
            % CHECK ANGLE OF ARRIVAL OF SEQUENCE
        elseif button=='a',
            try
                length_temp=0.1*ceil(abs(diff(current)*10));
                dtagaudit_synch_int_angle(tag,min(current),length_temp);
            catch
                disp(' An error occurred during dtagaudit_synch_int_aangle and action was aborted ')
            end
            
            % INSERT INSTANTANEOUS CUE AT CURSOR POSITION
        elseif button=='l',
            ss = input(' Enter comment... ','s') ;
            RES.cue = [RES.cue;[gx 0]] ;
            RES.stype{size(RES.cue,1)} = ss ;
            save tagaudit_RECOVER RES
            plotRES(AXc,RES,[tcue tcue+NS]) ;
            
            % DELETE OR CHANGE CUE AT CURSOR POSITION
        elseif button=='x',
            % Find all cues around synch cursor
            kres =(find(gx>=RES.cue(:,1)-0.1 & gx<sum(RES.cue')'+0.1)) ;
            if ~isempty(kres),
                % Find the type of cues and allow users to edit them
                ktype = RES.stype(kres) ;
                ktype=inputdlg(ktype,'Edit/remove cues',[1,40],ktype);
                % If cancel button is pressed, ktype will be empty cell
                if ~isempty(ktype)
                    % Otherwise, go through cell array of new sound types
                    for i=1:length(ktype)
                        RES.stype(kres(i))=ktype(i);
                        kempty(i)=isempty(ktype{i});
                    end
                    
                    % If there are empty sound types, remove them
                    kempty=kres(find(kempty));
                    kkeep = setxor(1:size(RES.cue,1),kempty) ;
                    RES.cue = RES.cue(kkeep,:) ;
                    RES.stype = {RES.stype{kkeep}} ;
                end
                plotRES(AXc,RES,[tcue tcue+NS]) ;
                clear kempty kres
                
            else
                fprintf(' No saved cue at cursor\n') ;
            end
            
            % GO FORWARD ONE SCREEN
        elseif button=='f',
            tcue = tcue+floor(NS)-0.5 ;
            done = 1 ;
            
            % GO BACK ONE SCREEN
        elseif button=='b',
            tcue = max([0 tcue-NS+0.5]) ;
            done = 1 ;
            
            % PLAY ALL SOUND
        elseif button=='p',
            if ~isempty(bs),
                xf = filter(bs,as,x) ;
                sound(volume*xf,afs/SOUND_DF,16) ;
            else
                sound(volume*x,afs/SOUND_DF,16) ;
            end
            
            % PLAY SHORT SELECTION AND SAVE TEMP FILE
        elseif button=='i',
            if all(current>0),
                length_temp=0.1*ceil(abs(diff(current)*10));
                [x_temp,afs_org] = dtagwavread(tag,min(current),length_temp);
                if isempty(x), return, end
                
                % First, query for sound label
                if QUERY_LABEL,
                    label=inputdlg({'Sound label'},'Add label (skA, skB, etc) to save',[1,40],{''});
                else
                    label=[];
                end
                if ~isempty(label),
                    outputname = ['call types/' tag '_' num2str(round(min(current))) '_' char(label) '.wav'];
                    %outputname = [tag '_' num2str(round(min(current))) '.wav'];
                else
                    outputname = 'tempsound.wav';
                end
                
                % Then try to save
                try
                    wavwrite(x_temp(:,CH),afs_org,16,outputname);
                catch
                    disp([' Error writing to ' outputname ', close other programs and try again'])
                end
                
                % Finally, filter and play sound
                if ~isempty(bs),
                    xf = filter(bs,as,x_temp) ;
                    xf = resample(xf,48e3,afs_org);
                    sound(volume*xf,48e3/SOUND_DF,16) ;
                else
                    xf = resample(x_temp,48e3,afs_org);
                    sound(volume*xf,48e3/SOUND_DF,16) ;
                end
            else
                fprintf('Invalid click: Need to mark interval before pressing i to play sound from this interval\n')
            end
            
            % ADJUST SYNCH AUDITING PARAMETERS
        elseif button=='o',
            % Find out which tags are possible
            if isempty(possibletags)
                letters='abcdefghij'; n=0;
                % Check for other tags from same day
                for i=1:length(letters),
                    % Rewrite this to find all possible tags, figure out when
                    % they were on and off the animal, and find estimated
                    % time difference - then keep track of this in structure
                    if tagver==2,
                        othertagfile = makefname([tag(1:8) letters(i)],'AUDIO',1);
                    elseif tagver==3,
                        othertagfile = d3makefname([tag(1:8) letters(i)],'AUDIO',1);
                    end
                    slashes = [union(findstr(othertagfile,'/'),findstr(othertagfile,'\'))];
                    othertagfile = othertagfile(1:slashes(end));
                    if exist(othertagfile) % compatible w d2
                        if ~strcmp(tag,[tag(1:8) letters(i)]) % Exclude focal tag
                            n=n+1;                            % Add tag to possible tag list
                            possibletags{n}=[tag(1:8) letters(i)];
                        end
                    end
                end
            end
            
            % Choose tags for comparing
            InitValue=[];
            if ~isempty(othertags)
                for i=1:length(othertags)
                    if strmatch(char(othertags(i)),possibletags,'exact')
                        InitValue=[InitValue strmatch(char(othertags(i)),possibletags,'exact')];
                    end
                end
            end
            [kd,ok] = listdlg('ListString',possibletags,'SelectionMode','multiple',...
                'InitialValue',InitValue,'OKString','Accept','Name','Tags',...
                'PromptString',['Select tags to compare with focal ' tag],...
                'ListSize',[300 100]);
            
            if ~ok,
                continue
            end
            othertags = possibletags(kd) ;
            
            % Check time shifts available for each tag
            if ~(length(othertags)==length(time_shift))
                time_shift=zeros(1,length(othertags));
            end
            
            % Adjust initialvalue
            clear InitValue
            for i=1:length(othertags)
                InitValue{i}=num2str(time_shift(i));
            end
            
            % Adjust time shifts for each tag
            kd=inputdlg(othertags,'Time offset',[1,30],InitValue);
            for i=1:length(othertags)
                if ~isnan(str2double(kd{i}))
                    time_shift(i)=str2double(kd{i});
                else
                    time_shift(i)=0;
                end
            end
            
            % Save synch audit settings for quick retrieval
            synch_fname = [tag(1:9) '_synch.mat'] ;
            save (synch_fname,'othertags','time_shift')
            
        elseif button=='r',
            if isempty(othertags)
                disp(' First adjust synch parameters using "o" ')
                continue
            end
            
            try
                length_temp=0.1*ceil(abs(diff(current)*10));
                dtagaudit_synch_int_revise(tag,min(current),othertags,length_temp,time_shift);
            catch
                disp(' An error occurred during tagaudit_synch_int_revise and action was aborted ')
                disp([' Trying to compare ' tag ' with the following tags:' ])
                disp(char(othertags))
            end
            
            
            % SYNCH SOUND COMPARISON WITH OTHER TAGS
        elseif button=='t'
            if isempty(othertags)
                disp(' First adjust synch parameters using "o" ')
                continue
            end
            
            if length(othertags)<6
                try
                    tic
                    [AXS_synch]=dtagaudit_synch_int(tag,tcue,RES,othertags,time_shift,NS,AFS_RESAMPLE);
                    toc
                catch
                    disp(' An error occurred during d3tagaudit_synch_int and action was aborted ')
                    disp([' Trying to compare ' tag ' with the following tags:' ])
                    disp(char(othertags))
                end
            else
                disp('d3tagaudit_synch currently supports comparing focal tag with up to 4 other tags')
            end
            
            
            % SYNCH ANALYSIS OF ANGLE-OF-ARRIVAL ON MULTIPLE TAGS
        elseif button=='T'
            if isempty(othertags)
                disp(' First adjust synch parameters using "o" ')
                continue
            end
            
            try
                dtagaudit_synch_int_aoa(tag,tcue,RES,othertags,NS,time_shift,AFS_RESAMPLE);
            catch
                disp(' An error occurred during d3tagaudit_synch_aoa_int and action was aborted ')
                disp([' Trying to compare ' tag ' with the following tags:' ])
                disp(char(othertags))
            end
            
            % MARK SEQUENCE ON NON-FOCAL TAGS
        elseif button=='S'
            if isempty(othertags)
                disp(' Cannot mark multiple sequences. First adjust synch parameters using "o" ')
                continue
            end
            
            % Check if simultaneous spectrogram plot and handles exist
            if isempty(findobj('type','figure','name','Simultaneous Spectrogram Plot'))
                disp(' Cannot mark multiple sequences. Need to display simultaneous spectrograms ')
                continue
            elseif isempty(AXS_synch)
                disp(' Cannot mark multiple sequences. Need to display simultaneous spectrograms first')
                continue
            end
            
            
            % Call that figure
            figure(findobj('type','figure','name','Simultaneous Spectrogram Plot'))
            
            % Set up artificial cursors in each synch plot, make them
            % invisible until changed to visible
            for i=1:length(AXS_synch(:,1))
                axes(AXS_synch(i,1));
                hold on, hhh_synch(i) = plot([0 0],0.8*min([afs/2000 MAXYONFREQDISPLAY])*[1 1],'k*-','visible','off') ; hold off
            end
            
            MARK=1;
            current_synch = [0 0] ;
            axis_selected = [0 0] ;
            
            % Activate relevant axis to bring out crosshairs from ginput
            axes(AXS_synch(1,1));
            
            while MARK
                
                % Find two points
                [gx_synch gy_synch button_synch] = ginput(1) ;
                
                % Verify that user clicked in right figure
                if gcf~=findobj('type','figure','name','Simultaneous Spectrogram Plot'),
                    continue,
                end
                
                if button_synch==1,
                    if gy_synch<0 | gx_synch<0 | isempty(find(AXS_synch(:,1)==gca))
                        disp('Invalid click')
                    else
                        % Update selection
                        current_synch=[current_synch(2) gx_synch];
                        axis_selected=[axis_selected(2) gca];
                        % Make all synch cursors invisible
                        set(hhh_synch,'visible','off')
                        % Then activate the cursor that is active and update position
                        thistag=find(AXS_synch(:,1)==axis_selected(2));
                        set(hhh_synch(thistag),'visible','on','XData',current_synch(find(axis_selected==gca)),'YData',0.8*min([afs/2000 MAXYONFREQDISPLAY])*ones(length(find(axis_selected==gca)),1));
                    end
                    
                    
                elseif button_synch=='s'
                    
                    if diff(axis_selected) | all(axis_selected==0)
                        disp('Click twice in same axis before marking sequence')
                    else
                        
                        % Find axis that was clicked in and corresponding tag
                        thistag=find(AXS_synch(:,1)==axis_selected(2));
                        
                        % Decide cue type
                        if thistag==1,
                            ss=char(inputdlg(tag,'Cue type',[1,40]));
                            if ~isempty(char(ss))
                                cc = sort([current_synch]) ;
                                RES.cue = [RES.cue;[cc(1) diff(cc)]] ;
                                RES.stype{size(RES.cue,1)} = ss ;
                                save tagaudit_RECOVER RES
                                plotRES(AXc,RES,[tcue tcue+NS]); % Update auditing plot
                                plotRES(AXS_synch(thistag,2),RES,[tcue tcue+NS]); % Update synch plot
                            end
                        else
                            % Load audits from other tags,
                            % add cue to cue list, and save again
                            ss=char(inputdlg(othertags(thistag-1),'Cue type',[1,40]));
                            if ~isempty(char(ss))
                                cc = sort([current_synch]) ;
                                RESTEMP=loadaudit(char(othertags(thistag-1)));
                                RESTEMP.cue = [RESTEMP.cue;[cc(1) diff(cc)]] ;
                                RESTEMP.stype{size(RESTEMP.cue,1)} = ss ;
                                saveaudit(char(othertags(thistag-1)),RESTEMP) ;
                                plotRES(AXS_synch(thistag,2),RESTEMP,get(AXS_synch(thistag,2),'XLIM')) ;
                                clear RESTEMP
                                fclose('all') % Force matlab to close any files that are still open
                            end
                        end
                    end
                    
                    
                elseif button_synch=='i'
                    
                    if diff(axis_selected) | all(axis_selected==0) | diff(current_synch)>NS
                        disp('Click twice in same axis before playing sequence')
                    else
                        
                        % Find axis that was clicked in and corresponding tag
                        thistag=find(AXS_synch(:,1)==axis_selected(2));
                        
                        % Define what tag recording is selected
                        if thistag==1,
                            play_tag = tag ;
                        else
                            play_tag = char(othertags(thistag-1)) ;
                        end
                        
                        if all(current_synch>0),
                            length_temp=0.1*ceil(abs(diff(current_synch)*10));
                            [x_temp,afs_org] = dtagwavread(play_tag,min(current_synch),length_temp) ;
                            if isempty(x), return, end
                            
                            % First, query for sound label
                            if QUERY_LABEL,
                                label=inputdlg({'Sound label'},'Add label (skA, skB, etc) to save',[1,40],{''});
                            else
                                label=[];
                            end
                            if ~isempty(label),
                                outputname = ['call types/' play_tag '_' num2str(round(min(current_synch))) '_' char(label) '.wav'];
                                %outputname = [play_tag '_' num2str(round(min(current_synch))) '.wav'];
                            else
                                outputname = 'tempsound.wav';
                            end
                            
                            % Then, write original sound extract to disk
                            try
                                wavwrite(x_temp(:,CH),afs_org,16,outputname);
                            catch
                                disp([ ' Error writing to ' outputname ', close other programs and try again'])
                            end
                            
                            % Finally, play back sound
                            if ~isempty(bs),
                                xf = filter(bs,as,x_temp) ;
                                xf = resample(xf,48e3,afs_org);
                                sound(volume*xf,48e3/SOUND_DF,16) ;
                            else
                                xf = resample(x_temp,48e3,afs_org);
                                sound(volume*xf,48e3/SOUND_DF,16) ;
                            end
                        else
                            fprintf('Invalid click: Need to mark interval before pressing i to play sound from this interval\n')
                        end
                        
                    end
                    
                elseif button_synch=='l'
                    
                    % Cancel if selection is not within a valid axis
                    if isempty(find(AXS_synch(:,1)==gca))
                        disp('Click in valid spectrogram axis to activate that axis first')
                        continue
                    end
                    
                    % Buttons other than left clicking does not update selected
                    % axis - instead, check that click is within selected axis
                    axislimits=get(axis_selected(2),'YLim');
                    if gy_synch<axislimits(1) | gy_synch>axislimits(2)
                        disp('To change audit cue, press x while hovering in active spectrogram')
                        continue
                    end
                    
                    % Find axis that was clicked in and corresponding tag
                    thistag=find(AXS_synch(:,1)==axis_selected(2));
                    
                    if thistag==1,
                        ss=char(inputdlg(tag,'Cue type',[1,40]));
                        if ~isempty(char(ss))
                            RES.cue = [RES.cue;[gx_synch 0]] ;
                            RES.stype{size(RES.cue,1)} = ss ;
                            save tagaudit_RECOVER RES
                            plotRES(AXc,RES,[tcue tcue+NS]) ;
                            plotRES(AXS_synch(thistag,2),RES,[tcue tcue+NS]); % Update synch plot
                        end
                    else
                        % Load audits from other tags,
                        % add cue to cue list, and save again
                        ss=char(inputdlg(othertags(thistag-1),'Cue type',[1,40]));
                        if ~isempty(char(ss))
                            RESTEMP=loadaudit(char(othertags(thistag-1)));
                            RESTEMP.cue = [RESTEMP.cue;[gx_synch 0]] ;
                            RESTEMP.stype{size(RESTEMP.cue,1)} = ss ;
                            saveaudit(char(othertags(thistag-1)),RESTEMP) ;
                            plotRES(AXS_synch(thistag,2),RESTEMP,get(AXS_synch(thistag,2),'XLIM')) ;
                            clear RESTEMP
                        end
                    end
                    
                    % Delete cue at cursor
                elseif button_synch=='x',
                    
                    % Buttons other than left clicking does not update selected
                    % axis - instead, check that click is within selected axis
                    axislimits=get(axis_selected(2),'YLim');
                    if gy_synch<axislimits(1) | gy_synch>axislimits(2)
                        disp('To change audit cue, press x while hovering in active spectrogram')
                        continue
                    end
                    
                    % Find axis that was clicked in and corresponding tag
                    thistag=find(AXS_synch(:,1)==gca);
                    
                    % If tag is focal, RES is current audit structure
                    if thistag==1,
                        % Find all cues around synch cursor
                        kres =(find(gx_synch>=RES.cue(:,1)-0.1 & gx_synch<sum(RES.cue')'+0.1)) ;
                        if ~isempty(kres),
                            % Find the type of cues and allow users to edit them
                            ktype = RES.stype(kres) ;
                            ktype=inputdlg(ktype,'Edit/remove cues',[1,40],ktype);
                            % If cancel button is pressed, ktype will be empty cell
                            if ~isempty(ktype)
                                % Otherwise, go through cell array of new sound types
                                for i=1:length(ktype)
                                    RES.stype(kres(i))=ktype(i);
                                    kempty(i)=isempty(ktype{i});
                                end
                                
                                % If there are empty sound types, remove them
                                kempty=kres(find(kempty));
                                kkeep = setxor(1:size(RES.cue,1),kempty) ;
                                RES.cue = RES.cue(kkeep,:) ;
                                RES.stype = {RES.stype{kkeep}} ;
                                
                                % Update information in single audit and
                                % simultaneous audit screen
                                plotRES(AXc,RES,[tcue tcue+NS]) ;
                                plotRES(AXS_synch(thistag,2),RES,[tcue tcue+NS]); % Update synch plot
                            end
                            clear kempty kres
                            
                        else
                            fprintf(' No saved cue at cursor\n') ;
                        end
                        
                        % If tag is non-focal, load audit as RESTEMP
                    else
                        % WARNING: BACK UP AUDIT TXT FILES AND TEST
                        % THOROUGHLY ANY CHANGES MADE TO THIS SECTION
                        
                        % Load audits from other tag,
                        RESTEMP=loadaudit(char(othertags(thistag-1)));
                        % Find all cues around cursor
                        kres =(find(gx_synch>=RESTEMP.cue(:,1)-0.1 & gx_synch<sum(RESTEMP.cue')'+0.1)) ;
                        if ~isempty(kres),
                            % If cancel button is pressed, ktype will be empty cell
                            ktype = RESTEMP.stype(kres) ;
                            ktype=inputdlg(ktype,'Edit cues (empty to remove)',[1,30],ktype);
                            if ~isempty(ktype)
                                for i=1:length(ktype)
                                    % Otherwise, go through cell array of new sound types
                                    RESTEMP.stype(kres(i))=ktype(i);
                                    kempty(i)=isempty(ktype{i});
                                end
                                % Remove empty sound types
                                kempty=kres(find(kempty));
                                kkeep = setxor(1:size(RESTEMP.cue,1),kempty) ;
                                RESTEMP.cue = RESTEMP.cue(kkeep,:) ;
                                RESTEMP.stype = {RESTEMP.stype{kkeep}} ;
                                % Now save audit again
                                saveaudit(char(othertags(thistag-1)),RESTEMP) ;
                                % Update synch plot
                                plotRES(AXS_synch(thistag,2),RESTEMP,get(AXS_synch(thistag,2),'XLIM'));
                            end
                            clear kempty kres
                        else
                            fprintf(' No saved cue at cursor\n') ;
                        end
                        clear RESTEMP
                    end
                    
                    % If other button is pressed, including q, go back to single auditing screen
                else
                    disp(' Aborting nonfocal call marking')
                    MARK=0; continue
                end
            end
            
            
            
            % ADJUST AUDIT SELECTION
        elseif button==1,
            if gy<0 | gx<tcue | gx>tcue+NS
                fprintf('Invalid click: commands are f b s l p x q\n')
                
            else
                current = [current(2) gx] ;
                set(hhh,'XData',current) ;
                if ~isempty(p),
                    fprintf(' -> %6.1f\t\tdiff to last = %6.2f\t\tp = %6.1f\t\tfreq. = %4.2f kHz\n', ...
                        gx,diff(current),p(round(gx*fs)),gy) ;
                else
                    fprintf(' -> %6.1f\t\tdiff to last = %6.1f\t\tfreq. = %4.2f kHz\n', ...
                        gx,diff(current),gy) ;
                end
            end
        end
    end
end



function   [cc,xx] = findclicks(x,thresh,fs,fh)
%
%     clicks = findclicks(x,thresh,fs,fh)
%     Return the time cue to each click in x, in seconds
%     x is a signal vector
%     thresh is the click detection threshold - try a value of 0.01 to
%     begin with
%     fs is the sampling rate of x [default = 32000]
%     fh is the high pass filter frequency [default = 10000]
%
%     mark johnson, WHOI
%     August, 2000
%     modified June 2003

% for sperm whales use:
tc = 2.5e-3 ;           % power averaging time constant in seconds

% for ziphius use:
tc = 0.5e-3 ;           % power averaging time constant in seconds

blanking = 20e-3 ;      % blanking time after a click is detected before another

[b a] = cheby1(6,0.5,fh/fs*2,'high') ;
pp = 1/tc/fs ;

xf = filter(b,a,x);
xx = filter(pp,[1 -(1-pp)],abs(xf)) ;
cc = [] ;

if thresh==0,
    return
end

cc = find(diff(xx>thresh)>0)/fs ;
done = 0 ;

if isempty(cc),
    return ;
end

while ~done,
    kg = find(diff(cc)>blanking) ;
    done = length(kg) == (length(cc)-1) ;
    cc = cc([1;kg+1]) ;
end
return


function plotRES(AXc,RES,XLIMS)

axes(AXc)
if ~isempty(RES.cue),
    kk = find(sum(RES.cue')'>XLIMS(1) & RES.cue(:,1)<=XLIMS(2)) ;
    if ~isempty(kk),
        plot([RES.cue(kk,1) sum(RES.cue(kk,:)')']',0.2*ones(2,length(kk)),'k*-') ;
        for k=kk',
            text(max([XLIMS(1) RES.cue(k,1)+0.1]),0.6,RES.stype{k},'FontSize',10) ;
        end
    else
        plot(0,0,'k*-') ;
    end
else
    plot(0,0,'k*-') ;
end

set(AXc,'XLim',XLIMS,'YLim',[0 1]) ;
bc = get(gcf,'Color') ;
set(AXc,'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;
return