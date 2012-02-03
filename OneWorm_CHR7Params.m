%% LENIENT VARIABLES TO BE SPECIFIED  - modifiy these to fit your application
vers='chronos7INTENSITY(8_9_11)- single cell flouresence'

%% ImageAnalysisMode
framerate=29
MicroM_Pixel=1

EvenImgBgSub= 'single'
imgfmt='*jpg'
%path smoothng options
SmoothInt=5% 
SmoothMeth='moving'
%Smoothstart=25 %dont plot until smoothed index is this long <<specify
%within Plot4StacktiffOneWorm < later import all params as single matirx
poshead=[]


SnglImgProofMd = 'off'  % _LEVE THIS OFF MOSTLY _single image proofing MODE 
DataCapMode = 'StackGiff' %'StackTiff','Simple'
allow_img= 'n'%allow_img= input ('Allow images? (Y/N)', 's')
proofingImgVIS ='off'%'off' Visualize scored stacks as they are built
ProofingImgs = 'n' % extra proofing images
stoppoint='n'


 App='New'
 FntSz = 14;
invertImage= 'y'

%% more parameters - DEACTIVATED OPTIONS
IsolatePlate = 'N' %for standard chronos Apps this will be 'Y', Need a plate in img1 and img2    
AlignImMod='n'
StringentFilter = 'n'%** STARTS OFF AS NOto run first filter the stringent filter will kick in if you get < 2.1* MaxWorms (2 points for each moving worm)
%>startdate=datenum('16-Jun-2011 12:00:00')
halveimage= 'n' %use to avoid before after worm fusion, and check for robutness of timepoiints
CropEdge = 'n' % turn this on if edges are giving you a hard time > there is a crop Manually

%% PARTICLE COUNT filters -Lenient variable set
MaxFilt=2  %number of filters tried before continuing
MaxWormFactor=1 %skips if MaxWorms * MaxWormFactor is exceeded 
MaxWorms=1 % MAXIMUM PARTICLES COUNTED
MinWorms=1
%% THRESHOLDING FILTERS value -OR use graythesh to dynamically determine     
dynamicTH = 'y' %when set to automatic seems to switch between .0001 (bad) and .49 (ok)
thresh_hold=.49; % dynamicTH deactivates this one .25
dynamicBndLim = 'prc' %'stdv', 'prc''static' 
%>>NumStd=7 %three standard deviations do not seem to do it (7 matches a .25% boundlimit (about))
BndLim=.65;  % PRC MODE- dynamicBndLim deactivates this one .25
%>>Minbnd=-20% STATIC MODE- specifically remove the middle 'bnd' % of colors
%>>Maxbnd=20 

%% PARTICLE FILTERS
FiltApp='SZ_Ax_BB'%'all'
LowLim=500; UpLim=1200; % AREA %col10 <<NEEDED TO 
MajAxL=50; MajAxU=300; %col 8   <<NEEDED TO RAise 
MinAxL=5; MinAxU=300; %col 9   <<NEEDED TO RAise <<needed to lower Min 10>5
BndBx4L=25; BndBx4U=1000; %col14 <<NEEDED TO RAise 
TotAxL = 85; TotAxU = 200  %(MajAxL+MinAxL)

%% Second pass- STRINGENT variable set
if strcmpi('y', StringentFilter) 

    dynamicTH = 'y'
    thresh_hold=.49; % dynamicTH deactivates this one

    %BW threshold value -OR use graythesh to dynamically determine
   %> MaxElapseTm=14 %maximum (hrs) interval to allow subtraction (prevents physical artafacts)
   %> MinTime_elapsed= 0.01 %minimum time I cant see movement with the short time intervals

    dynamicBndLim = 'n'
    MaxElapseTm=15 %maximum (hrs) interval to allow subtraction (prevents physical artafacts)
    BndLim=.75;  %***<< maybe .35?
    LowLim=11.5; UpLim=35; 
    eccL=.7; eccU=.99; %col 5 
    MajAxL=4; MajAxU=15;  %<<<lowered for partial worms
    MinAxL=1.5; MinAxU=5; 
    ExtL=.25; ExtU=.8; % optional parameter input, OTHERWISE USE DEFAULTS
    BndBx4L=2.9; BndBx4U=12; 
end    

%% get working foler if required

%WorkFldTest = exist ('DataDir') % if this is running without wormviewsuite you need to specify directory
%if WorkFldTest > 0 
%else
%message= ('choose the "working" folder') % this is the folder that will
%collect the results
%DataDir = uigetdir('choose the "working" folder')
%end

%imgSize = 'L' %>>imgSize= input ('which image size? (Small (S), Large (L)', 's') 
%>Programs = ('/Applications/MATLAB_R2010a.app/toolbox/ObjectIntensityModule')


%% DEFAULT filter set 
%radSqr=780000; %masks a circle around area of the agar
if exist ('BndLim') == 0; BndLim=.35; end %specify +/- values around average to be removed.
if exist ('LowLim') == 0; LowLim=15; end %<<<LowLim=18.9, UpLim=45 
if exist ('UpLim') == 0; UpLim=250; end
if exist ('eccL') == 0; eccL=.2; end
if exist ('eccU') == 0; eccU=1; end %<eccL=.7, eccU=1;   %eccentricity filter e.g. linear objects: 
if exist ('MajAxL') == 0; MajAxL=10; end
if exist ('MajAxU') == 0; MajAxU=200; end
if exist ('MinAxL') == 0; MinAxL=2; end
if exist ('MinAxU') == 0; MinAxU=20; end
if exist ('ExtL') == 0; ExtL=0; end
if exist ('ExtU') == 0; ExtU=1; end %ranges from 0-1
if exist ('BndBx4U') == 0; BndBx4U=20; end
if exist ('BndBx4L') == 0; BndBx4L=7; end 
            
 