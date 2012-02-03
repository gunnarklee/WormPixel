function GetWorm(varargin)

%OneWorm_CHRONOS7.m
%GetWorm - start point for One Worm analysis pipeline

% G. Kleemann (8/3/11)
%from twoimageDropCHR7_3Robot.m

%Version History G. Kleemann

%take OneWorm and size of an object from a stack of images.
%modified from 2 image drop CHRONOS_Robot
%Set up for PNG files now, and folders starting with 'pics_'

%twoImagDrop pulled from Chronos structure and rewritten as a standalone
%app (7/6/09) -GAK

%twoImagDrop modified to add non aligning and cropping version for Jarlaths
%protocol where backgroud is generated then subtraced from a single image (3/21/10) -GAK

%Robot version 7_2 optimized for cox robot
%select ellipse mask.

%new features 7/2011 - to reduce failures nad allow troubleshooting
%1-separately score before and after to avoud particle loss by lumping 20
%ositions into one worm _also one may be more reliable

%2-upper amd lower time difference limits for pairwise subtraction
%Need xceptions for early time points since contamination will be low and
%worms moving well <not yet implimented>
%Need forward search for pair with mintime diff <not yet implimented>


%3-LOW QUALITY IMAGE PAIR FILTERS
%maximum worm Number error
%filters with stringent filters switchrd in and pairs
%dropped when the count is very inflated (seems to happpen once in a while
%with othetr pairs stil well scored

%Particle detction error 
%drop the image if whole plate is the largest particle

% Parser for case when inputdir and outputdir are specified
p = inputParser;
p.addRequired('inputdir', @isdir);
p.addRequired('outputdir', @isdir);
p.addOptional('trialname', 'trialname', @ischar);

% Parse Inputs
try
    p.parse(varargin{:})
catch e1
    try
        inputdir = uigetdir('choose the "Source" folder');
        outputdir = uigetdir('choose the "Destination" folder');
        p.parse(inputdir, outputdir, varargin{:})
    catch e2
        exception = MException(...
            'twoimageDropCHR7_3Robot:arglist',...
            'Error in input argument list');
        exception = addCause(exception, e1);
        exception = addCause(exception, e2);
        throw(exception)
    end
end
DataDir = p.Results.inputdir;
Alldata = p.Results.outputdir;
TrialName = p.Results.trialname;

disp(['Input directory: ' DataDir]);
disp(['Output directory: ' Alldata]);
disp(['Trial name: ' TrialName]);

%% Make new directories
RUNfinalDir = [Alldata, 'RESULTS', filesep, TrialName, 'RUNfinal'];
ErrorDir = [Alldata, 'RESULTS', filesep, TrialName, 'ErrorDir'];
mkdir(RUNfinalDir);
mkdir(ErrorDir)
ProcessDate=(date);
%>>>save ([RUNfinalDir, filesep, 'Params'], 'Params');
%% SPECIFY FILTER PARAMETERS 

%OPTIONAL Filter Param Specification OR revert to defaults
%thresh_hold=.0001; %BW threshold value -OR use graythesh to dynamically determine
%BndLim=.25; LowLim=8; UpLim=35; eccL=.2; eccU=1; MajAxL=5; MajAxU=200; MinAxL=0; MinAxU=20; ExtL=0; ExtU=1; % optional parameter input, OTHERWISE USE DEFAULTS

%cd([CodeDir, '/Chronos7/ObjectOneWormModule'])
OneWorm_CHR7Params  % load parameters an prefs
%twoImgDpPramStruct;

%% Setup variables
%% set up variables
scrsz = get(0,'ScreenSize');
AliveData={};
ImPrProcTm=[];
centroidLs=[];
HeadPosnLs=[];
Skiplog={};
OneWorm_Data={};
Image_PropertiesAll={};
SpineStack=[]
CurveMtx=[]


%% GET FINAL.mat file names if they exist
%>>[DateFldrNms RecentFldr namedate] = GetFolderMat (DataDir, 'final', 'final')

%% specify directories
% Alldata = DataDir;
%     %cd (Alldata)
%     OutDir = [Alldata '/RESULTS'];
%     OutDir2=[Alldata filesep TrialName 'RUNfinal'];  
%     mkdir (OutDir); mkdir(OutDir2);
% %set default values _should work if image format is correct
%     %dtPtls = 'dot';
    %dtPtMark=dtPtls; 
    
%% GET folder Names
dirOutput = dir(fullfile(Alldata, 'PIC_*')); %specifiy surce folder
    DateFldrNms = {dirOutput.name}';
%% check each folder or crop params if absent get them  
crop_position = [75 10 425 425];

for W=1:length(DateFldrNms)
    
    dirOutputCropPar = dir(fullfile([Alldata, '/',DateFldrNms{W}], 'CropParam.mat')); %list the images
    if size( dirOutputCropPar,1) < 1 % if there are no crop parameters, then get them
    DateFldrNms2 = {dirOutputCropPar.name}';
    
    dirOutputtif = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
    if size(dirOutputtif,1) < 1
        dirOutputtif = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
    end
    
    %cd([Alldata, '/',DateFldrNms{W}]);
    imgtmp=imread([Alldata filesep DateFldrNms{W} filesep dirOutputtif(5).name]); %load the first image in the folder
    imgtmp=rgb2gray(imgtmp);
    clear ('posctr'); clear ('posedge');clear('xctr');            %imgtmp = img1; imsubLC = EvnImgSubt;
            %imgtmp2 = img2; imsubLC2 = img2;
%>>            cd([CodeDir, '/Chronos7/Chronos_2ImageDrop']);  
            %CropCircleDropCHR7; %imgCP1=imgCP; 
             %cd(CodeDir); 
             CropRectOneWorm; % collects elipse object in posctr
            cd([Alldata, '/',DateFldrNms{W}]);
            %save ('CropParam.mat', 'radSqr', 'xctr', 'mask', 'center'); %save parameters
            save ([Alldata filesep DateFldrNms{W} filesep], 'CropParam.mat', 'mask', 'posctr'); %save parm -posctr s ellipse
        close all %clean up 
    else 
        %cd([Alldata, '/',DateFldrNms{W}]);
        load([Alldata filesep DateFldrNms{W} filesep dirOutputCropPar.name]);
    end
end  

%%
%Single image mode
if strcmpi('off', SnglImgProofMd);
FldMax=length(DateFldrNms); FldStart=1;
else %single image on
FldStart=FldMax;   
end

%%
for W=FldStart:FldMax; %loop single images
      
    centroidLs=[];
    HeadPosnLs=[];
    BBratio=[];
    centrSmthY=[];
   
    dirOutput2 = dir(fullfile([Alldata filesep DateFldrNms{W}], imgfmt)); %list the images
    if size(dirOutput2 ,1) < 1
        dirOutput2  = dir(fullfile([Alldata filesep DateFldrNms{W}], imgfmt)); %list the images
    end
%  sort the structure into dates order from time stamp    
        %convert dates to serial numbers
        for n= 1:size(dirOutput2,1);
        dirOutput2(n).datenum=datenum({dirOutput2(n).date});
        end
        
        [unused, order] = sort([dirOutput2(:).datenum]); %sort by number
        dirOutput2 = dirOutput2(order);
        {dirOutput2.date}';

    DateFldrNms2 = {dirOutput2.name}';
    %cd([Alldata, '/',DateFldrNms{W}]);


close all
Images=[];

    if strcmpi(EvenImgBgSub, 'y') % only for non-simple mode   
        CumulImg=double(imread(DateFldrNms2{1}));
    end    
%%
%read the images into the work space as a STACK
%in EvenImgBgSub mode, make cumulimages.
for ImN=1:size(DateFldrNms2,1);
if strcmpi('off', SnglImgProofMd); % switch for single pairmode
    ImgNmMax=size(DateFldrNms2,1)-1;
    ImgStrt=1;
else
    ImgStrt=ImgNmMax;
end    
    
   load([Alldata filesep DateFldrNms{W} filesep 'CropParam.mat']);
   startimageName=[DateFldrNms2{1}];
   imageName=[DateFldrNms2{ImN}];
   suffix=DateFldrNms{W}(findstr('_',DateFldrNms{W})+1:end);
   startImgCt=str2num(startimageName(findstr(suffix, startimageName)+length(suffix):findstr('.', startimageName)-1));
   CurrimageCt=str2num(imageName(findstr(suffix, imageName)+length(suffix):findstr('.', imageName)-1));
   imageCt=CurrimageCt-startImgCt+1;
  
    timeintv =imageCt/framerate;

    img1=imread([Alldata filesep DateFldrNms{W} filesep DateFldrNms2{ImN}]);
    img1=rgb2gray(img1);
    subimg=double(img1);

close all 

%MASK OFF THE PERIMITER
img=(subimg.*mask); 
lng1=length(img(1,:));
lng2=length(img(:,1));

%% limit search neighborhood to
%LAST WORM BOUNDING BOX
if ImN > 2

%Search only in a padded region boinding the last worm
PadPrc=1.5
%PaddedBox = [boundingBox(1)*PadPrc,boundingBox(2).*PadPrc,boundingBox(3)*PadPrc*4,boundingBox(4)*PadPrc*4]   
%PlotBoundBox(img, PaddedBox) % to check..
PlotBoundBox(img, boundingBox) % to check..
%PaddedBox,boundingBox
BBmask = zeros(size(img));
%Be sure that padded BB does not stretch the mask beyond the edges of img 
BBmask (boundingBox(2)-boundingBox(4).*PadPrc:boundingBox(2)+boundingBox(4).*2.*PadPrc,...
        boundingBox(1)-boundingBox(3).*PadPrc:boundingBox(1)+boundingBox(3).*2.*PadPrc)= ones;
%Be sure that padded BB does not stretch the mask beyond the edges of img 
BBmask=BBmask(1:size(img,1), 1:size(img,2))  

img=(img.*BBmask); %% pass the masked image into the particle analysis to avoit off target particles. 
figure; imshow(uint8(img));
%stoppt=input('next step?', 's') 
%clear ('PaddedBox', 'boundingBox')
end 

% MASK BY INTensity ** removing middle tones
%be careful this can lead to rescaling and washing
% out of extreeme features. this seems to start with BndLim < .8...
% Mask=(img < (bnd+avg) | img > (avg-bnd)); %threshold set as bounds around AVERAGE
%.9 masks almost nothing, .00001 masks everything
% Mask=(img > (bnd+avg) | img < (avg-bnd)); %this masks the extreeme values around bndlim
% OR .9 masks almost everything and, .00001 masks almost nothing
close all
%BndLim=.3

% dynamically determine bound limits from image characteristics
    switch dynamicBndLim; 
        case 'stdv'    
        Minbnd=-NumStd*(std2(img));  %3 std is not nearly enough
        Maxbnd=NumStd*(std2(img)); 
        case 'prc'
        Minbnd=BndLim*(min(min(img))); %specifically remove the middle 'bnd' % of colors
        Maxbnd=BndLim*(max(max(img)));
        case 'static'
        Minbnd=Minbnd; %specifically remove the middle 'bnd' % of colors
        Maxbnd=Maxbnd;  
        
    end 

%avg=mean2(img); doesnt seem to work in matlab 2010
avg=mean(img(:));
    %>>>[xm,ym]=meshgrid(1:lng1,1:lng2); %meshgrid adapts to picture size
    Mask=(img < (Minbnd) | img > (Maxbnd)); %threshold set as bounds around AVERAGE
    %>>Mask=(img > (bnd+avg) | img < (avg-bnd)); %this masks the extreeme
    %values around bndlim
    
        if strcmpi(EvenImgBgSub, 'y');% only for non-simple mode 
        Mask=imcomplement(Mask); % off for straight subtract ON for adjusted subtract
        end
    % get STATS for the mask
    imgRange=[(max(max(img))); min(min(img)); range(range(img));avg; BndLim];
    RangeRmvd=[Minbnd;Maxbnd];
    totPixels= size(Mask, 1) *size(Mask, 2);
    %Mask=(img > (avg+bnd)| img < (avg-bnd)); % perhaps substiute 0 for (avg-bnd) here?
    pxlmsk=totPixels-sum(sum(Mask));
    PrcMasked=(pxlmsk./totPixels).*100;
      %Mask=Mask-1; % masked values come in as ones so change them to zeros
            if (strcmpi (allow_img, 'y')); 
                figure;imagesc(Mask); title ('Mask image');
            end
      Mask=double(Mask);
      img=double(img);
      MasImg=(img.*Mask);
            if (strcmpi (allow_img, 'y'))       
%>            figure; imagesc(Straightsubimg); title ('Straightsubimg')
            figure; imagesc(MasImg); title ('MasImg')
            figure; imagesc(Mask); title ('Mask')
            figure ('position', scrsz); subplot (2,2,1); imagesc(subimg); title ('original image'); colorbar; subplot (2,2,2); imagesc(Mask);title ('Mask image'); colorbar; subplot (2,2,3); imagesc(MasImg); title ('MasImg image'); colorbar;
            end   


%% IDENTIFY OBJECTS and FILTER DATA
%%BW THRESHLDING - the boundlimit subtraction above takes care of the thresholding      
            
                   imgBW=abs(double(MasImg));
                   if strcmpi(invertImage, 'y'); imgBW=imcomplement(imgBW); end
                   % thresh_hold=.0001 - thresholding from .0001 -.999 is
                   % making no differnece!
                        if (strcmpi(dynamicTH, 'y'));
                        thresh_hold = graythresh(imgBW);% dynamically determine threshold
                        end
                        imgBW=im2bw(imgBW,thresh_hold);% apply threshold

                    if strcmpi(allow_img, 'y'); figure; imagesc(imgBW); end
            %BW=imgBW; 
            %GetImgPropsDropAug11;
            clear ('Image_PropertiesAll', 'F', 'AlabeldAll', 'gnumAll', 'yy', 'yyy', 'xx', 'xxx', 'ym', 'xm', 'mask')
            [imgBWL, F, Image_PropertiesAll] = GetImgProps (imgBW, allow_img);
            
            if size(Image_PropertiesAll, 1) < 1 %values absent - make one dummy line
            Image_PropertiesAll= ones(1,21);
            end
            
%if strcmpi(allow_img, 'y')
%figure;
%subplot (2,2,2); imagesc(BW); title ('thresholded image'); 
%subplot (2,2,3); imshow(BW);colorbar;
%subplot (2,2,4); imagesc(imgBWL); title ('BW labeled image');%axis equal;
%subplot (2,2,1); imagesc(MasImg); title ('masked image');%axis equal;
%end

%% Filter and present data. >>** expand this ti include separate filtering
%% and tabulaiton for before/ afterlists
            imgdsp=img1; %J1_exp; 
            imgdsp2=imgBWL; 
            if strcmpi ('y', allow_img); figure ; imshow(imgBWL); end 
            clear ('Img_Propfilt')
            
         AppFltsAug11 %needs to add fliter on bounding box BndBxFlt4L.* BndBxFlt4U.*
         FltrandDisp4OneWorm
%% COUNT ERROR-  Too many particles, then SKIP THE PAIR
                if strcmpi('n', countgood)
                    %adjust your particle thresholds']); break
                    save ([ErrorDir filesep imageName(1:end-4), 'CountError.mat'], 'F',... %'Imagesfilt', 
                   'filtVal','imgBWL', 'img', 'img1', 'Img_Propfilt', 'Image_PropertiesAll');%'ProcessDate'
                    continue
                end
                
%% spine worm
            [SpineData, poshead2] = SpineWorm(Imagesfilt, Img_Propfilt, img1, ErrorDir, allow_img, stoppoint, poshead)    
            
            [poshead]=updatePoshead (poshead2, poshead)
            
            %check for errors        
              if strcmpi('n', SpineData.spinegood)
                    save ([ErrorDir filesep imageName(1:end-4), 'SpineError.mat'], 'F',... %'Imagesfilt', 
                   'filtVal','imgBWL', 'img1', 'Img_Propfilt', 'Image_PropertiesAll');%'ProcessDate'
                    continue
              end 
 
            FltNm2={};
            FltNm2= ['sze'; 'ecc'; 'mjx'; 'mnx'; 'ext'];
            FltPram2= [LowLim UpLim; eccL eccU; MajAxL MajAxU; MinAxL MinAxU; ExtL ExtU];

            % list of image attributes
            if strcmpi ('y', allow_img); figure; imagesc(imgBWL); end
            
            numObj=length (Img_Propfilt (:,1));

%some variables are repeated in the output
%there are two centroid pairs and two areas but dont remove them since the
%program refereces specific columns throughout

%these are graphed in "GetImgPropsDrop.m"
%for a description of a specific measure look up "regionprops" in help
Allsizes = Img_Propfilt(:,10);
sizesScored = Image_PropertiesAll(:,10); 

ImgPropHeaders =...
{'Centroid'; 'Centroid'; 'Area'; 'objectnumber'; 'Eccentricity';...
'Centroid'; 'Centroid'; 'MajorAxisLength'; 'MinAxisLength'; 'Area';... 
'BoundingBox'; 'BoundingBox'; 'BoundingBox'; 'BoundingBox'; 'Extent';... 
'EquivalentDiameter'; 'EulerNumber'; 'Perimeter'; 'perim/area'; 'majAx/minAx';...
'(majAx/minAx)/extent'};

    % Output work to .matfile put these in a single large floder with the run
    % date indicated 
      %cd(Alldata);
          
        
          %save (['final.mat'], 'ProcessDate', 'filtVal', 'FltNm2' ,'imgBWL',...
          %'imgCP',  'Img_Propfilt', 'Image_PropertiesAll','imsubLC', 'img1', 'img2', 'radSqr')
          
          save ([RUNfinalDir filesep imageName(1:end-4), 'final.mat'],  'Imagesfilt', 'ImgPropHeaders', 'F',...
          'filtVal', 'FltNm2' ,'imgBWL', 'img1', 'Img_Propfilt', 'Image_PropertiesAll', 'SpineData');%'ProcessDate'
      
        %  save ([RUNfinalDir filesep imageName(1:end-4), 'final.mat'], 'img1','Image_PropertiesAll',...
        % 'Image_PropertiesList', 'Images', 'Imagesfilt', 'Img_Propfilt') %, 'CumulImg'      

majAx_minAx = (Img_Propfilt(1,20));
majA_minAx_extent = (Img_Propfilt(1,21));
CentrX = (Img_Propfilt(1,1));
CentrY = (Img_Propfilt(1,2));
areacell = (Img_Propfilt(1,3));



OneWorm_Data = [OneWorm_Data;  {DateFldrNms{W}, ceil(numObj),...
    areacell, areaboundingBox, majAx_minAx, majA_minAx_extent,... 
    CentrX, CentrY}];


              %AliveData={AbsTimeDays, timeintv, DateFldrNms{W}, numObj}
              %cellwrite(DataDir/, AliveData, 2)
%          
end % IMAGES LOOP

OneWorm_DataHeaders =  ['filename',',', 'ObjecCount',',',...             
'areacell',',','areaboundingBox',',','majAx_minAx',',', 'majA_minAx_extent',... 
    ',','CentrX',',', 'CentrY'];


save ([RUNfinalDir, filesep ,DateFldrNms{W},'_',DateFldrNms2{ImN}, 'AllData.mat']) % OutDir
save([RUNfinalDir, filesep ,DateFldrNms{W} 'OneWormdata.mat'], 'OneWorm_DataHeaders', 'OneWorm_Data')
end % FOLDERS LOOP
end

%% make sure you have a good head position
function [poshead]=updatePoshead (poshead2, poshead)
            if isempty(poshead2)
                poshead=poshead
            else 
                poshead=poshead2
            end 
end

