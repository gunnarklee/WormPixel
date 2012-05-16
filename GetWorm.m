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

ProcessDate=(date);
%% SPECIFY FILTER PARAMETERS

%OPTIONAL Filter Param Specification OR revert to defaults
%thresh_hold=.0001; %BW threshold value -OR use graythesh to dynamically determine
%BndLim=.25; LowLim=8; UpLim=35; eccL=.2; eccU=1; MajAxL=5; MajAxU=200; MinAxL=0; MinAxU=20; ExtL=0; ExtU=1; % optional parameter input, OTHERWISE USE DEFAULTS

OneWorm_CHR7Params  % load parameters an prefs


%% Setup variables
scrsz = get(0,'ScreenSize');
AliveData={};
ImPrProcTm=[];
centroidLs=[];
HeadPosnLs=[];
Skiplog={};
OneWorm_Data={};
Image_PropertiesAll={};
SpineStack=[];
CurveMtx=[];



%% GET folder Names
dirOutput = dir(fullfile(Alldata, 'PIC_*')); %specifiy source folder
DateFldrNms = {dirOutput.name}';
if isempty(DateFldrNms); error('No "PIC_..." folder, Check folder name and paths'); end

%% for each "PIC" folder check to see that the required files are there
% If not already done Crop the useful area from the picture
CheckGetCrop(Alldata, DateFldrNms, imgfmt)

% If filter params are missing, have user identify 5 worms and get params
Particleparams (Alldata, DateFldrNms, imgfmt, dynamicTH, thresh_hold)


%% get head position

if strcmpi('off', SnglImgProofMd); %Single image mode
    FldMax=length(DateFldrNms); FldStart=1;
else %single image on
    FldStart=FldMax;
end

for W=FldStart:FldMax; %loop folders
    %% new directory for each video
    %% Make new directories
    RUNfinalDir = [Alldata, 'RESULTS', filesep, TrialName, filesep, DateFldrNms{W} 'RUNfinal'];
    ErrorDir = [Alldata, 'RESULTS', filesep, TrialName, filesep, DateFldrNms{W} 'ErrorDir'];
    mkdir(RUNfinalDir);
    mkdir(ErrorDir);
    
    %New Partilce parameters
    load([Alldata filesep DateFldrNms{W} filesep 'FltrParams.mat'])
    
    %update filter values
    LowLim=FltrParams.ParticleFilt.LowLim;
    UpLim=FltrParams.ParticleFilt.UpLim;
    MajAxL=FltrParams.ParticleFilt.MajAxL;
    MajAxU=FltrParams.ParticleFilt.MajAxU;
    MinAxL=FltrParams.ParticleFilt.MinAxL;
    MinAxU=FltrParams.ParticleFilt.MinAxU+(MinAxU*.8);  %<<<raised to avoud dropping
    TotAxU=FltrParams.TotAxU;
    TotAxL=FltrParams.TotAxL;
    
    centroidLs=[];
    HeadPosnLs=[];
    BBratio=[];
    centrSmthY=[];
    
    dirOutput2 = dir(fullfile([Alldata filesep DateFldrNms{W}], imgfmt)); %list the images
    if size(dirOutput2 ,1) < 1
        error('I can not find any images')
    end
    
    
    %% sort the structure into dates order from time stamp
    %if strcmpi(SortByName, 'y')
    %    [dirOutput2] = SortByName (dirOutput2)
    % sometimes this mis-orders files if the time stamp is wrong
    %end
    
    DateFldrNms2 = {dirOutput2.name}'; % cell array to matrix
    
    close all
    Images=[];
    poshead=[];
    for ImN=1:size(DateFldrNms2,1);
        
        varStruct=[] ;
        CropPar=load([Alldata filesep DateFldrNms{W} filesep 'CropParam.mat']);%reload each time
        mask =CropPar.mask; posctr=CropPar.posctr;
        if strcmpi('off', SnglImgProofMd); % switch for single pairmode
            ImgNmMax=size(DateFldrNms2,1)-1;
            ImgStrt=1;
        else
            ImgStrt=ImgNmMax;
        end
        
        startimageName=[DateFldrNms2{1}];
        imageName=[DateFldrNms2{ImN}];
        suffix=DateFldrNms{W}(findstr('_',DateFldrNms{W})+1:end);
        startImgCt=str2num(startimageName(findstr(suffix, startimageName)+length(suffix):findstr('.', startimageName)-1));
        CurrimageCt=str2num(imageName(findstr(suffix, imageName)+length(suffix):findstr('.', imageName)-1));
        imageCt=CurrimageCt-startImgCt+1;
        
        timeintv =imageCt/framerate;
        
        img1=imread([Alldata filesep DateFldrNms{W} filesep DateFldrNms2{ImN}]);
        img1=imresize(img1, 2);
        if isrgb(img1); img1=rgb2gray(img1); end
        
        %% Remove blotchy background
        
        if strcmpi (smoothbkg, 'y')
            SE=strel('disk',10);
            imgTH=imbothat(img1, SE);
            imgTH=imcomplement(imgTH);
            imgAdj=adapthisteq(imgTH); %Opt?
            if (strcmpi (allow_img, 'y'))
                figure; imshow(imgTH)
                figure; imshow(img1)
                figure; imshow(imgAdj)
            end
            subimg=double(imgAdj);
        else
            subimg=double(img1);
        end
        
        %%
        close all
        
        %MASK OFF THE PERIMITER
        img=(subimg.*mask);
        lng1=length(img(1,:));
        lng2=length(img(:,1));
        
        %% limit search neighborhood to
        %LAST WORM BOUNDING BOX
        if ImN > 2
            %Search only in a padded region boinding the last worm
            PadPrc=1.5;
            
            %if there is no bounding box then make one
            %boundingBox=Img_Propfilt(:,11:14)
            if exist('boundingBox')==0; boundingBox=[200.  250.   100.0000   100.0000]; end
            
            BBmask = zeros(size(img));
            %Be sure that padded BB does not stretch the mask beyond the edges of img
            BBmask (boundingBox(2)-boundingBox(4).*PadPrc:boundingBox(2)+boundingBox(4).*2.*PadPrc,...
                boundingBox(1)-boundingBox(3).*PadPrc:boundingBox(1)+boundingBox(3).*2.*PadPrc)= ones;
            %Be sure that padded BB does not stretch the mask beyond the edges of img
            BBmask=BBmask(1:size(img,1), 1:size(img,2));
            
            img=(img.*BBmask); %% pass the masked image into the particle analysis to avoit off target particles.
            if (strcmpi (allow_img, 'y'))
                figure; imshow(uint8(img));
                PlotBoundBox(img, boundingBox); % to check..
            end
            %stoppt=input('next step?', 's')
            %clear ('PaddedBox', 'boundingBox')
        end
        
        %% MASK BY INTensity ** removing middle tones
        %be careful this can lead to rescaling and washing
        close all
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
        avg=mean(img(:));
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
            figure; imagesc(MasImg); title ('MasImg');
            figure; imagesc(Mask); title ('Mask');
            figure ('position', scrsz); subplot (2,2,1); imagesc(subimg); title ('original image'); colorbar; subplot (2,2,2); imagesc(Mask);title ('Mask image'); colorbar; subplot (2,2,3); imagesc(MasImg); title ('MasImg image'); colorbar;
        end
        
        close all
        
        %% IDENTIFY OBJECTS and FILTER DATA
        %%BW THRESHLDING - the boundlimit subtraction above takes care of
        %%the thresholding
        if strcmpi(invertImage, 'y'); MasImg=imcomplement(MasImg); end
        % thresh_hold=.0001 - thresholding from .0001 -.999 is
        % making no differnece!
        
        if (strcmpi(dynamicTH, 'y'));
            imgBW=abs(double(MasImg));
            thresh_hold = graythresh(imgBW);% dynamically determine threshold
            imgBW=imcomplement(im2bw(imgBW,thresh_hold));% apply threshold
            %figure; imagesc(imgBW);
        elseif (strcmpi(dynamicTH, 'y-ShortCircuit'));
            imgBW=abs(imgTH); %SKIPS OVER INTENSITY AND DOUBLE opreations
            thresh_hold = graythresh(imgBW);% dynamically determine threshold
            imgBW=imcomplement(im2bw(imgBW,thresh_hold));% apply threshold
            %figure; imagesc(imgBW);
        else
            imgBW=abs(double(MasImg));
            imgBW=im2bw(img1,thresh_hold);% apply threshold
            %figure; imagesc(imgBW);
        end
        
        if strcmpi(allow_img, 'y'); figure; imagesc(imgBW); end
        
        Image_PropertiesAll=[]; F=[]; AlabeldAll=[]; gnumAll=[]; yy=[]; yyy=[]; xx=[]; xxx=[]; ym=[];xm=[];mask=[];
        [imgBWL, F, Image_PropertiesAll] = GetImgProps (imgBW, allow_img);
        
        if size(Image_PropertiesAll, 1) < 1 %values absent - make one dummy line
            Image_PropertiesAll= ones(1,21);
        end
        
        %% Filter and present data. >>** expand this ti include separate filtering
        %% and tabulaiton for before/ afterlists
        imgdsp=img1; %J1_exp;
        imgdsp2=imgBWL;
        if strcmpi ('y', allow_img); figure ; imshow(imgBWL); end
        Img_Propfilt=[];
        
        
        %% AppFlts
        %needs to add fliter on bounding box BndBxFlt4L.* BndBxFlt4U.*
        
        %set up filters to identify worms - logical matricies of ones and zeros...
        szFltL= double(Image_PropertiesAll(:,10)>LowLim); %make size filter
        szFltU= double(Image_PropertiesAll (:,10)<UpLim);
        EccFltL= double(Image_PropertiesAll (:,5)>eccL); %rows less than 1280 are in chanel 1
        EccFltU= double(Image_PropertiesAll (:,5)<eccU);% the eccentricity of the spot should be less than .8 (usually ~.5)
        MajAxFltL= double(Image_PropertiesAll (:,8)>MajAxL); %rows less than 1280 are in chanel 1
        MajAxFltU= double(Image_PropertiesAll (:,8)<MajAxU);% the eccentricity of the spot should be less than .8 (usually ~.5)
        MinAxFltL= double(Image_PropertiesAll (:,9)>MinAxL); %rows less than 1280 are in chanel 1
        MinAxFltU= double(Image_PropertiesAll (:,9)<MinAxU);% the eccentricity of the spot should be less than .8 (usually ~.5)
        ExtFltL= double(Image_PropertiesAll (:,15)>ExtL); %rows less than 1280 are in chanel 1
        ExtFltU= double(Image_PropertiesAll (:,15)<ExtU);
        TotAxFltU = (Image_PropertiesAll (:,8))+(Image_PropertiesAll (:,9))<TotAxU;
        TotAxFltL = (Image_PropertiesAll (:,8))+(Image_PropertiesAll (:,9))>TotAxL;
        
        
        
        
        
        FltrandDisp4OneWorm
        %% COUNT ERROR-  Too many particles, then SKIP THE PAIR
        if strcmpi('n', countgood)
            varStruct.filters.filtVal=filtVal;
            varStruct.images.img=img;
            varStruct.images.img1=img1;
            varStruct.analysis.F=F;
            varStruct.analysis.Image_PropertiesAll=Image_PropertiesAll;
            varStruct.filters.szFltL= szFltL;
            varStruct.filters.szFltU= szFltU;
            varStruct.filters.MajAxFltL= MajAxFltL; %rows less than 1280 are in chanel 1
            varStruct.filters.MajAxFltU= MajAxFltU;% the eccentricity of the spot should be less than .8 (usually ~.5)
            varStruct.filters.MinAxFltL= MinAxFltL; %rows less than 1280 are in chanel 1
            varStruct.filters.MinAxFltU=MinAxFltU;
            
            saveThis([ErrorDir filesep imageName(1:end-4), 'CountError.mat'], varStruct)
            continue
        end
        
        %% get head position for first worm
        
        [WmImgPad]=GetPadImg (pad, Imagesfilt) ;
        %WmImgPad=imresize(logical(WmImgPad), 2);
        
        if isempty(poshead)
            %display a few images to tell which part is the head
            Flipbook([Alldata filesep DateFldrNms{W}], DateFldrNms2(1:10));
            %function [poshead, WmImgPad]=GetHeadPosPad (pad, Imagesfilt,)
            mssg='drag point to head and double click';
            [poshead]=GetHead (poshead, WmImgPad, mssg);
            varStruct.Pos.poshead=poshead;
            close all
        end
        %% spine worm
        figure; subplot(1,2,1); imshow(img1); subplot(1,2,2);  imshow(WmImgPad);
        
        [SpineData, poshead2] = SpineWorm (WmImgPad, img1, allow_img, stoppoint, poshead, numpts);
        %
        [poshead]=updatePoshead (poshead2, poshead);
        close all
        
        %% BADSPINE DUMP
        %check for errors
        if strcmpi('n', SpineData.spinegood)
            varStruct.filters.filtVal=filtVal;
            varStruct.images.img=img;
            varStruct.images.img1=img1;
            varStruct.analysis.F=F;
            varStruct.analysis.Image_PropertiesAll=Image_PropertiesAll;
            varStruct.analysis.Img_Propfilt=Img_Propfilt;
            
            saveThis ([ErrorDir filesep imageName(1:end-4), 'SpineError.mat'],varStruct);
            
            continue
        end
        
        FltNm2={};
        FltNm2= ['sze'; 'ecc'; 'mjx'; 'mnx'; 'ext'];
        FltPram2= [LowLim UpLim; eccL eccU; MajAxL MajAxU; MinAxL MinAxU; ExtL ExtU];
        
        % list of image attributes
        if strcmpi ('y', allow_img); figure; imagesc(imgBWL); end
        
        numObj=length (Img_Propfilt (:,1));
        
        %these are graphed in "GetImgPropsDrop.m"
        %for a description of a specific measure look up "regionprops" in
        %help
        Allsizes = Img_Propfilt(:,3);
        sizesScored = Image_PropertiesAll(:,3);
        
        ImgPropHeaders =...
            {'Centroid'; 'Centroid'; 'Area'; 'objectnumber'; 'Eccentricity';...
            'Centroid'; 'Centroid'; 'MajorAxisLength'; 'MinAxisLength'; 'Area';...
            'BoundingBox'; 'BoundingBox'; 'BoundingBox'; 'BoundingBox'; 'Extent';...
            'EquivalentDiameter'; 'EulerNumber'; 'Perimeter'; 'perim/area'; 'majAx/minAx';...
            '(majAx/minAx)/extent'};
        
        % Output work to .matfile put these in a single large floder with the run
        % date indicated
        %cd(Alldata);
        
        %% PREPARE for separate save function (required for parallel processing)
        
        varStruct.filters.filtVal=filtVal;
        varStruct.filters.FltNm2= FltNm2;
        varStruct.images.imgBWL=imgBWL;
        varStruct.images.img=img;
        varStruct.images.img1=img1;
        varStruct.images.Imagesfilt=Imagesfilt;
        varStruct.analysis.F=F;
        varStruct.analysis.ImgPropHeaders=ImgPropHeaders;
        varStruct.analysis.Image_PropertiesAll=Image_PropertiesAll;
        varStruct.analysis.Img_Propfilt=Img_Propfilt;
        varStruct.SpineData=SpineData;
        
        SaveImNm=imageName(1:end-4);
        saveThis([RUNfinalDir filesep SaveImNm, 'final.mat'], varStruct);%'ProcessDate'
        
        
        majAx_minAx = (Img_Propfilt(1,20));
        majA_minAx_extent = (Img_Propfilt(1,21));
        CentrX = (Img_Propfilt(1,1));
        CentrY = (Img_Propfilt(1,2));
        areacell = (Img_Propfilt(1,3));
        boundingBox=Img_Propfilt(:,11:14);
        
        
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
end % FOLDERS LOO P
end %end main function

%% make sure you have a good head position
function [poshead]=updatePoshead (poshead2, poshead)
if isempty(poshead2)
    poshead=poshead;
else
    poshead=poshead2;
end
end

function CheckGetCrop(Alldata, DateFldrNms, imgfmt)
%% check each folder or crop params if absent get them
crop_position = [75 10 425 425];
for W=1:length(DateFldrNms)
    dirOutputCropPar = dir(fullfile([Alldata filesep DateFldrNms{W}], 'CropParam.mat')); %list the images
    if size( dirOutputCropPar,1) < 1 % if there are no crop parameters, then get them
        DateFldrNms2 = {dirOutputCropPar.name}';
        dirOutputtif = dir(fullfile([Alldata filesep DateFldrNms{W}], imgfmt)); %list the images
        if size(dirOutputtif,1) < 1
            dirOutputtif = dir(fullfile([Alldata filesep DateFldrNms{W}], imgfmt)); %list the images
        end
        %cd([Alldata, '/',DateFldrNms{W}]);
        imgtmp=imread([Alldata filesep DateFldrNms{W} filesep dirOutputtif(5).name]); %load the first image in the folder
        imgtmp=imresize(imgtmp, 2);
        if isrgb(imgtmp)
            imgtmp=rgb2gray(imgtmp);
        end
        CropRectOneWorm; % collects elipse object in posctr
        save ([Alldata filesep DateFldrNms{W} filesep 'CropParam.mat'], 'mask', 'posctr'); %save parm -posctr s ellipse
        close all %clean up
    end
end

end

%% get Param range for the worm
function Particleparams (Alldata, DateFldrNms, imgfmt, dynamicTH, thresh_hold)
%% get particle filter values by asking the user to idenfity the worm in 5 pictures.
WormProps=[];

for W=1:length(DateFldrNms) % for each folder check for filter params
    dirOutputCropPar = dir(fullfile([Alldata, '/',DateFldrNms{W}], 'FltrParams.mat')); %list the images
    if size( dirOutputCropPar,1) < 1 % if there are no filter params, then get them
        %>DateFldrNms2 = {dirOutputCropPar.name}'; % get image list
        dirOutputtif = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
        
        %Get 5 images at equal intervals across iamge set
        
       QueryImgs=[1:ceil(size(dirOutputtif,1)/6):ceil(size(dirOutputtif,1))]
       for imgnm=QueryImgs %load 5 images to choose particles
            imgtmp=imread([Alldata filesep DateFldrNms{W} filesep dirOutputtif(imgnm).name]); %load the first image in the folder
            imgtmp=imresize(imgtmp, 2);
            if isrgb(imgtmp); imgtmp=rgb2gray(imgtmp); end
            
            if strcmpi(dynamicTH, 'y')
                thresh_hold = graythresh(imgtmp);% dynamically determine threshold
            end
            
       %% Get thresholded image
            % allow threshold adjustment in first image
            if imgnm==1
                TH = 'y';
                while strcmpi (TH, 'n') == 0
                    imgBW=imcomplement(im2bw(imgtmp,thresh_hold));% apply threshold
                    figure; imagesc(imgBW);
                    display(thresh_hold)
                    TH=input('rethreshold? - value or "n"','s');
                    if strcmpi(TH, 'n')==0
                        thresh_hold=str2num(TH);%=input('new value')
                        dynamicTH = 'n';
                    end
                end
                display(thresh_hold)
            else 
                imgBW=imcomplement(im2bw(imgtmp,thresh_hold));% apply threshold  
            end
            
            %find particles and get partilce properties
            [imgBWL, F, Image_PropertiesAll] = GetImgProps (imgBW, 'y');
            if imgnm==1; StartPos=[ceil(size(imgBWL(:,:,1), 2)*.5), ceil(size(imgBWL (:,:,1), 1)*.5)]; end
            
            %identify the worm
            mssg='Double click on the worm centorid'
            [pos] = GetPoint(imgBWL, StartPos, mssg);
            StartPos=pos;
            %get the particle with the closest centroid to the selected location
            [row, mindiff]=CloseCentr(pos, Image_PropertiesAll);
            
            %capture list of worm props
            Currhit=Image_PropertiesAll(row,:);
            WormProps=[WormProps;Currhit];
        end
        
        %build and save new parameter spec sheet
        %% BUILD THE NEW PARTICLE FILTERS
        [upper, lower, stats]= GetParamLimits(WormProps(:,3), .4);
        FltrParams.ParticleFilt.LowLim=lower;
        FltrParams.ParticleFilt.UpLim=upper;  %col10 <<NEEDED TO RAise AREA TO <85 for clump; <5 for smallest only
        
        %MAJOR and MINOR AXIS NEED TO BE MORE PERMISSIVE
        [upper, lower, stats]= GetParamLimits(WormProps(:,8), .8);
        FltrParams.ParticleFilt.MajAxL=lower;
        FltrParams.ParticleFilt.MajAxU=upper; %col 8   <<NEEDED TO RAise AREA TO >24 for clump;  <3.3 for smallest only
        
        [upper, lower, stats]= GetParamLimits(WormProps(:,9), .8);
        FltrParams.ParticleFilt.MinAxL=lower;
        FltrParams.ParticleFilt.MinAxU=upper; %col 9   <<NEEDED TO RAise AREA TO >12 for clump;
        
        FltrParams.TotAxL=FltrParams.ParticleFilt.MinAxL+FltrParams.ParticleFilt.MajAxL;
        FltrParams.TotAxU=FltrParams.ParticleFilt.MinAxU+FltrParams.ParticleFilt.MajAxU;
        %% filters
        FiltApp='SZ_Ax_BB';%'all'
        
        save ([Alldata filesep DateFldrNms{W} filesep 'FltrParams.mat'], 'FltrParams', 'pos'); %save parm -posctr s ellipse
        close all %clean up
    end
end %folder loop
end

function [row, mindiff]=CloseCentr(pos, Image_PropertiesAll)

distancediff=[];
for Ptnm=1:size(Image_PropertiesAll, 1)
    
    %get each distance between the seleced position and each centroid
    currdist=pdist2(Image_PropertiesAll(Ptnm,1:2),pos);
    [row,col]=find(distancediff==min(distancediff));
    distancediff=[distancediff;currdist];
end
%row and dist for the closest particle
[row,col]=find(distancediff==min(distancediff));
mindiff=distancediff(row,:);
end

% get parameter limits from a list of a property
function [upper, lower, stats] = GetParamLimits(DataList, PRCbracket)
stats.mean=mean(DataList);
stats.range=mean(DataList);
stats.upper=stats.mean+stats.mean*PRCbracket;upper=stats.upper;
stats.lower=stats.mean-stats.mean*PRCbracket;lower=stats.lower;
stats.std=std(DataList);

end


function [dirOutput2] = SortByName (dirOutput2)

% sometimes this mis-orders files if the time stamp is wrong
% Make it optional

%convert dates to serial numbers >
for n= 1:size(dirOutput2,1);
    dirOutput2(n).datenum=datenum({dirOutput2(n).date});
end
%sort by date serial number the datetamp is wrong on some of these?
%'N2 A1 1_frame_0200.jpg', ... 0400.jpg, 0700.jpg, 1000.jpg but not 0300.jpg
%may have smting to do with when the image was saved?
[~, order] = sort({dirOutput2(:).name});
dirOutput2 = dirOutput2(order);
end


function [poshead]=GetHead (poshead, WmImgPad, mssg)
if isempty(poshead)
    [poshead] = GetPoint(WmImgPad, [ceil(size(WmImgPad (:,:,1), 2)*.85), ceil(size(WmImgPad (:,:,1), 1)*.85)], mssg);
end
end


function [WmImgPad]=GetPadImg (pad, Imagesfilt)
if iscell(Imagesfilt)
    WmImg=Imagesfilt{1,1};
else
    WmImg=Imagesfilt;
end
%pad image
[WmImgPad] = padImg (WmImg, pad);
end

%%>>>In Progress<<<
%>>function [Struct]=MakeStruct(varargin)

%>>p = inputParser;
%p.addRequired('inputdir', @isdir);
%p.addRequired('outputdir', @isdir);
%p.addOptional('trialname', 'trialname', @ischar);

% Parse Inputs
%>>p.parse(varargin{:})  %%%%"Must Be Scalar



%>>filtVal
%>>for n=1:length(elelmentList)
%>>    Element='F'
%>>    item=F

%>>Struct.(Element)=item



%>>end
%>>end


