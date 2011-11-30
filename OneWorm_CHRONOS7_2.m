%OneWorm_CHRONOS7.m
% G. Kleemann (8/3/11)

%take OneWorm and size of an object from a stack of images.
%modified from 2 image drop CHRONOS_Robot
%Set up for PNG files now, and folders starting with 'pics_'

%LongeMgrOveral3-
%3_15_08 manager file for Chronos1.2 by Gunnar Kleemann
%added opional start from cropped images (user option)
%added directories for Gunnars mac or lab PC (user option)
% in version 1.25 added dtPtMark which adjusts for different image namr
% markers.

%in chronos 2 (5/24/08) has been built to avoid frequent crashes due to
%memory errors.
%I have resructured the execution of the program so processes are executed in steps, only 
%parameters used for the next step are saved
% also, image pairs are in a single large folder to avoud recursve loops.


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

%Non-monotonic trend filter <not yet implimented>

%% GET FILE/SCRIPT INFORMATION 


cd (matlabroot)
cd toolbox
[stat,mess]=fileattrib('OneWormModule'); %get wormviews current attributes
%[stat,mess]=fileattrib('Chronos_2ImageDrop') %get two image drop's current attributes
MnProg = mess.Name; %dynamically specifiy the current womview locaiton
CodeDir = MnProg;
addpath(genpath(MnProg))
cd (MnProg)
scrsz = get(0,'ScreenSize');
AliveData={};
ImPrProcTm=[];
centroidLs=[];
HeadPosnLs=[];
Skiplog={};


message= ('choose the "Source" folder');
DataDir = uigetdir('choose the "Source" folder');

message= ('choose the "Destination" folder');
Alldata = uigetdir('choose the "Destination" folder');

%% SPECIFY FILTER PARAMETERS 

%OPTIONAL Filter Param Specification OR revert to defaults
%thresh_hold=.0001; %BW threshold value -OR use graythesh to dynamically determine
%BndLim=.25; LowLim=8; UpLim=35; eccL=.2; eccU=1; MajAxL=5; MajAxU=200; MinAxL=0; MinAxU=20; ExtL=0; ExtU=1; % optional parameter input, OTHERWISE USE DEFAULTS

    OneWorm_Data={};
    Image_PropertiesAll={};
  
    
TrialName=input('specify Trial Name e.g. "AIMtest10_10_2010"', 's');
%cd([CodeDir, '/Chronos7/ObjectOneWormModule'])
OneWorm_CHR7Params  % load parameters an prefs

%% specify directories
Alldata = DataDir;
    cd (Alldata)
    OutDir = [Alldata '/RESULTS'];
    OutDir2=[Alldata filesep TrialName 'RUNfinal'];  
    mkdir (OutDir); mkdir(OutDir2);
%set default values _should work if image format is correct
    %dtPtls = 'dot';
    %dtPtMark=dtPtls; 
    
%% GET folder Names
dirOutput = dir(fullfile(Alldata, 'PIC_*')); %specifiy surce folder
    DateFldrNms = {dirOutput.name}';
%% check each folder or crop params if absent get them   
for W=1:length(DateFldrNms)
    
    dirOutputCropPar = dir(fullfile([Alldata, '/',DateFldrNms{W}], 'CropParam.mat')); %list the images
    if size( dirOutputCropPar,1) < 1 % if there are no crop parameters, then get them
    DateFldrNms2 = {dirOutputCropPar.name}';
    
    dirOutputtif = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
    if size(dirOutputtif,1) < 1
        dirOutputtif = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
    end
    
    cd([Alldata, '/',DateFldrNms{W}]);
    imgtmp=imread(dirOutputtif(5).name); %load the first image in the folder
    imgtmp=rgb2gray(imgtmp);
    clear ('posctr'); clear ('posedge');clear('xctr');            %imgtmp = img1; imsubLC = EvnImgSubt;
            %imgtmp2 = img2; imsubLC2 = img2;
%>>            cd([CodeDir, '/Chronos7/Chronos_2ImageDrop']);  
            %CropCircleDropCHR7; %imgCP1=imgCP; 
             cd(CodeDir); CropRectOneWorm; % collects elipse object in posctr
            cd([Alldata, '/',DateFldrNms{W}]);
            %save ('CropParam.mat', 'radSqr', 'xctr', 'mask', 'center'); %save parameters
            save ('CropParam.mat', 'mask', 'posctr'); %save parm -posctr s ellipse
        close all %clean up 
    else 
        cd([Alldata, '/',DateFldrNms{W}]);
        load(dirOutputCropPar.name);
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
for W=FldStart:FldMax; %loop Well folders
      
    centroidLs=[];
    HeadPosnLs=[];
    BBratio=[];
    centrSmthY=[];
   
    dirOutput2 = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
    if size(dirOutput2 ,1) < 1
        dirOutput2  = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
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
    cd([Alldata, '/',DateFldrNms{W}]);


close all
Images=[];

    if strcmpi(EvenImgBgSub, 'y') % only for non-simple mode   
        CumulImg=double(imread(DateFldrNms2{1}));
    end    
%%
%read the images into the work space as a STACK
%in EvenImgBgSub mode, make cumulimages.
for ImN=1:size(DateFldrNms2,1);
ImN;
if strcmpi('off', SnglImgProofMd); % switch for single pairmode
    ImgNmMax=size(DateFldrNms2,1)-1;
    ImgStrt=1;
else
    ImgStrt=ImgNmMax;
end    
    
%%     
%>for ImN=ImgStrt:ImgNmMax % loop images
    %>dirOutputCP = dir(fullfile([Alldata, '/',DateFldrNms{W}], 'CropParam.mat'));
    cd([Alldata, '/',DateFldrNms{W}]);
    load('CropParam.mat');
   
    %check that the immage are not too separate temporally 
%%%%DOES NOT WORK FOR IMAGE STACKS YET SINCE IMAGES ARE GENERATED AFTER THE RUN
   %timeElapse = datevec(dirOutputtif(ImN).datenum -  dirOutputtif(1).datenum)
   %timeintv = (timeElapse(1,2)*2592000 + timeElapse(1,3)*86400 + timeElapse(1,4)*3600+ timeElapse(1,5)*60 + timeElapse(1,6))
   startimageName=[DateFldrNms2{1}];
   imageName=[DateFldrNms2{ImN}];
   suffix=DateFldrNms{W}(findstr('_',DateFldrNms{W})+1:end);
   %imageCt=str2num(imageName(length(imageName)-7:end-4))
   startImgCt=str2num(startimageName(findstr(suffix, startimageName)+length(suffix):findstr('.', startimageName)-1));
   CurrimageCt=str2num(imageName(findstr(suffix, imageName)+length(suffix):findstr('.', imageName)-1));
   imageCt=CurrimageCt-startImgCt+1;
   
   
   timeintv =imageCt/framerate;
   %timeintv = (timeElapse(1,2)*2592000 + timeElapse(1,3)*86400 + timeElapse(1,4)*3600+ timeElapse(1,5)*60 + timeElapse(1,6))
%
 %>   AbsTimeElapse=datevec(dirOutput2(ImN+1).datenum - startdate)
 %>   AbsTimeDays=(AbsTimeElapse(1,5)/3600 + AbsTimeElapse(1,4)/24 + AbsTimeElapse(1,3))

    img1=imread(DateFldrNms2{ImN});
    img1=rgb2gray(img1);
    subimg=double(img1);
    
cd (Alldata); dir % list the data so you can see format
%dtPtMark=dtPtls; % two varables used differntly

% Apply circular CROP then apply OneWorm mask.
close all 

%MASK OFF THE PERIMITER
img=(subimg.*mask); 
lng1=length(img(1,:));
lng2=length(img(:,1));

% MASK BY INTESNITY ** removing middle tones
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
    [xm,ym]=meshgrid(1:lng1,1:lng2); %meshgrid adapts to picture size
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


%%IDENTIFY OBJECTS and FILTER DATA

%%BW THRESHLDING - the boundlimit subtraction above takes care of the thresholding      
            
                   imgBW=abs(double(MasImg));
                   if strcmpi(invertImage, 'y'); imgBW=imcomplement(imgBW); end
                   % thresh_hold=.0001 - thresholding from .0001 -.999 is
                   % making no differnece!
                   % Dynamic thesshold option
                        if (strcmpi(dynamicTH, 'y'));
                        thresh_hold = graythresh(imgBW);% dynamically determine threshold
                        end
                        imgBW=im2bw(imgBW,thresh_hold);% apply threshold

                    if strcmpi(allow_img, 'y'); figure; imagesc(imgBW); end
            %>> clear ('Image_PropertiesAll')
            
        %%>>    cd([CodeDir, '/Chronos7/ObjectOneWormModule'])
            
            BW=imgBW; GetImgPropsDropAug11;
            
            %%if there are no values make one dummy line
            if size(Image_PropertiesAll, 1) < 1
            Image_PropertiesAll= ones(1,21);
            end
            
if strcmpi(allow_img, 'y')
%>figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]); 
figure;
subplot (2,2,2); imagesc(BW); title ('thresholded image'); 
subplot (2,2,3); imshow(BW);colorbar;
subplot (2,2,4); imagesc(imgBWL); title ('BW labeled image');%axis equal;
subplot (2,2,1); imagesc(MasImg); title ('masked image');%axis equal;
%if size (Imagesfilt, 1) > 1 
%figure;subplot(2,2,1);imshow(Imagesfilt{1});
%subplot(2,2,2);imshow(Imagesfilt{2}); title('I see more than one objects');
%end
end


%input ('is this OK')            
            
%% Filter and present data. >>** expand this ti include separate filtering
%% and tabulaiton for before/ afterlists
            imgdsp=img1; %J1_exp; 
            imgdsp2=imgBWL; 
            if strcmpi ('y', allow_img); figure ; imshow(imgBWL); end 
            clear ('Img_Propfilt')
            
         cd(CodeDir)
         AppFltsAug11 %needs to add fliter on bounding box BndBxFlt4L.* BndBxFlt4U.*
         FltrandDisp4OneWorm
%% COUNT FILTER Too many particle, then SKIP THE PAIR
                if strcmpi('n', countgood)
                    continue
                    %>>>>error(['detected more than one partlcle in image  ', num2str(ImN) ,'  adjust your particle thresholds']); break
                end  
            FltNm2={};
            FltNm2= ['sze'; 'ecc'; 'mjx'; 'mnx'; 'ext'];
            FltPram2= [LowLim UpLim; eccL eccU; MajAxL MajAxU; MinAxL MinAxU; ExtL ExtU];

            % list of image attributes
            if strcmpi ('y', allow_img); figure; imagesc(imgBWL); end
            
            %>>numObj=numObjOne+numObjTwo
            numObj=length (Img_Propfilt (:,1));
%>end
%          %>>% MasFltExtData
%
% RESULTs are in Img_Propfilt and Image_PropertiesAll
%Headers are in "ImgPropHeaders" 

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
      cd(Alldata);
          
        
          %save (['final.mat'], 'ProcessDate', 'filtVal', 'FltNm2' ,'imgBWL',...
          %'imgCP',  'Img_Propfilt', 'Image_PropertiesAll','imsubLC', 'img1', 'img2', 'radSqr')
          
          save ([OutDir2  filesep DateFldrNms{W},'_',DateFldrNms2{ImN}, 'final.mat'],  'Imagesfilt', 'ImgPropHeaders', 'F',...
          'filtVal', 'FltNm2' ,'imgBWL', 'img1', 'Img_Propfilt', 'Image_PropertiesAll');%'ProcessDate'
      

majAx_minAx = (Img_Propfilt(1,20));
majA_minAx_extent = (Img_Propfilt(1,21));
CentrX = (Img_Propfilt(1,1));
CentrY = (Img_Propfilt(1,2));
areacell = (Img_Propfilt(1,3));

OneWorm_Data = [OneWorm_Data;  {DateFldrNms{W}, ceil(numObj),...
    areacell, areaboundingBox, majAx_minAx, majA_minAx_extent,... 
    CentrX, CentrY}]


              %AliveData={AbsTimeDays, timeintv, DateFldrNms{W}, numObj}
              %cellwrite(DataDir/, AliveData, 2)
%          
end % IMAGES LOOP

OneWorm_DataHeaders =  ['filename',',', 'ObjecCount',',',...             
'areacell',',','areaboundingBox',',','majAx_minAx',',', 'majA_minAx_extent',... 
    ',','CentrX',',', 'CentrY'];


save ([OutDir, filesep ,DateFldrNms{W},'_',DateFldrNms2{ImN}, 'AllData.mat']) % OutDir
save([OutDir, filesep ,DateFldrNms{W} 'OneWormdata.mat'], 'OneWorm_DataHeaders', 'OneWorm_Data')
end % FOLDERS LOOP

