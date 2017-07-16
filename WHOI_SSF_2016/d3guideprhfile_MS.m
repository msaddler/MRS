% Generate DTAG3 PRH and CAL file
% Online GEOMAG converter: http://www.ngdc.noaa.gov/geomag-web/#igrfwmm
% MS 2016.06.20 EDITS --> Custom for blue whales, Chile

cd 'G:\Mark\code'
% SETUP PATHS % --> run startup.m from G:\tag
% Add path so MatLab can find Tag Tools
addpath ('G:\Mark\dtagtools\d2matlab')
addpath ('G:\Mark\dtagtools\d3analysis')
addpath ('G:\Mark\dtagtools\d3matlab_2014')
addpath ('G:\Mark\dtagtools\x3toolbox')
addpath ('G:\Mark\dtagtools\x3toolbox\XML4MATv2')
addpath ('G:\Mark\dtagtools\d3cal')
% Set paths for DTAG files
settagpath('cal', 'G:\Mark\data_bm\cal');
settagpath('raw', 'G:\Mark\data_bm\raw');
settagpath('prh', 'G:\Mark\data_bm\prh');
settagpath('audit', 'G:\Mark\data_bm\audit') ;
settagpath('audio', 'G:\Mark\data_bm\audio');

clear, close all

% Define expedition specific information stored in cal files
% author     = 'F. H. Jensen';
% email      = 'frants.jensen@gmail.com; plt@st-andrews.ac.uk';
% expedition = 'Sarasota 2016' ;
% DECL       = -5.28;   % magnetic declination in degrees, - = West, + = East
% INCL       = 56.58;   % magnetic inclination in degrees
% FSTRENGTH  = 45226; % total field strength, microTesla
% UTC2LOC    = -4 ;

% Define tag deployment constants
tag = 'bm14_076b' ;
recdir = ['G:\' tag(1:4) '\' tag] ; % this will create errors if drive is mapped to different drive
prefix = [tag(1:2) tag(6:end)] ;
df     = 1 ; % MUST BE SET: df=1

% % Define extra deployment specific information stored in cal file
% location       = 'NE of White Key';
% ANIMAL.ID      = 'FB209' ;
% ANIMAL.name    = 'C552' ;
% ANIMAL.gender  = 'Female' ;
% ANIMAL.description = '12-year old mother' ;
% ANIMAL.age     =  12;
% TAGON.TIME     = [2016 05 06 16 27 29] ; % Get this from xml field cue time
% TAGON.POSITION = [27.40395 82.64573] ;
% TAGON.PLACEMENT= 'Midway between blowhole and dorsal fin';

% FIND AND IMPLEMENT TAG CALIBRATION
[CAL,DEPLOY] = d3deployment ( recdir , prefix , tag ) ;

% READ SENSOR DATA
X   = d3readswv ( recdir , prefix , df ) ;

% Fix timing errors in sensor data
%X = d3fixsensorgap(X,df,DEPLOY);

% CONSIDER RUNNING fixsensormatrix.m

% PERFORM CALIBRATION OF TEMPERATURE AND PRESSURE
%etime(TAGON.RELEASE,TAGON.TIME)/3600
[p,CAL]=d3calpressure(X,CAL,'full') ;

% PERFORM CALIBRATION OF ACCELEROMETERS AND MAGNETOMETERS
%[A,CAL,fs] = d3calacc(X,CAL,'full',5); % <-- MS: If error occurs here, run in segments below!
[A,CAL,fs] = d3calacc(X,CAL,'bias',5);
[A,CAL,fs] = d3calacc(X,CAL,'p',5);
[A,CAL,fs] = d3calacc(X,CAL,'sens',5);
[A,CAL,fs] = d3calacc(X,CAL,'cross',5);

[M,CAL] = d3calmag(X,CAL,'full',5); % <-- Possibly run in segments also!


any(isnan(A))
any(isnan(M))
for j=1:3, A(1:120,j)=A(121,j);  M(1:120,j)=M(121,j); end, p(1:120)=p(121);
any(isnan(A))
any(isnan(M))


% check the lengths of the A, M and p matrices
S = [size(A,1),size(M,1),length(p)] ;
if length(unique(S))>1,
   n = min(S) ;
   fprintf(' Length mismatch in A, M and p. Trimming by %d samples\n',max(S)-n) ;
   A = A(1:n,:) ;
   M = M(1:n,:) ;
   p = p(1:n) ;
end


% SAVE PRH AND CAL DATA FOR AUDITING AND OTHER PROCESSING
saveprh(tag,'p','M','A','fs')
d3savecal(tag,'CAL',CAL);
disp('prh and cal files saved!')


% % ADD METADATA TO CAL FILE
% 
% % Who made files and recordings?
% d3savecal(tag,'AUTHOR',author);
% d3savecal(tag,'EMAIL',email);
% d3savecal(tag,'EXPEDITION',expedition);
% 
% % What were the magnetic field characteristics of the expedition site?
% d3savecal(tag,'DECL',DECL);
% d3savecal(tag,'INCL',INCL);
% d3savecal(tag,'FSTRENGTH',FSTRENGTH);
% d3savecal(tag,'UTC2LOC',UTC2LOC)
% 
% % What animal was tagged?
% d3savecal(tag,'ANIMAL',ANIMAL);
% d3savecal(tag,'PHOTO','N/A')
% 
% % Where and when was it tagged?
% d3savecal(tag,'LOCATION',location);
% d3savecal(tag,'TAGON',TAGON)
% 
% 
% % PERFORM TAG TO WHALE FRAME ROTATION OF ACC AND MAG
% % AND ESTIMATE PITCH/ROLL/HEADING
% 
% % Ideal orientation OTAB = [1 0 0 0 0]. 
% % Correct to actual OTAB found using prhpredictor
% PRH = prhpredictor (p,A,fs,4,2,'descent') ;
% OTAB = [1 0 0 0 0]; % ONLY IDEAL, i.e. NO TAG2WHALE ROTATION
% OTAB = [1 0 mean(PRH(:,2)) mean(PRH(:,3)) mean(PRH(:,4))]; %avg orientation from prhpredictor
% OTAB = [1 0 mean(PRH(:,2)) mean(PRH(:,3)) 0];
% [Aw,Mw] = tag2whale (A,M,OTAB,fs) ;
% d3savecal(tag,'OTAB',OTAB);
% 
% 
% % d3makeprhfile(recdir,prefix,df,tag); THIS DOES NOT WORK WITH MODIFIED X 
% % INSTEAD, DO MANUAL ROTATION AS FOLLOWS:
% 
% % report on trustworthiness of heading estimate
% dp = (A.*M*[1;1;1])./norm2(M) ;
% incl = asin(mean(dp)) ;     % do the mean before the asin to avoid problems
%                            % when the specific acceleration is large
% sincl = asin(std(dp)) ;
% fprintf(' Mean Magnetic Field Inclination: %4.2f\260 (%4.2f\260 RMS)\n',...
%        180/pi*incl, 180/pi*sincl) ;
% 
% % Calculate pitch, roll, heading, then correct for declination
% [pitch roll] = a2pr(Aw) ;
% [head vm incl] = m2h(Mw,pitch,roll) ;
% if exist('DECL'),
%    head = head + DECL*pi/180 ;      % adjust heading for declination angle in radians
% else
%     disp('Warning: No declination specified')
% end
% 
% % Save complete PRH file
% saveprh(tag,'p','A','M','fs','Aw','Mw','pitch','roll','head') ;
%    
% 
% R=loadaudit(tag);
% figure,plott(p,fs)

