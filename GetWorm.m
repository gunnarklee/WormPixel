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

%% Setup variables / matricies
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

ImgPropHeaders =...
    {'Centroid'; 'Centroid'; 'Area'; 'objectnumber'; 'Eccentricity';...
    'Centroid'; 'Centroid'; 'MajorAxisLength'; 'MinAxisLength'; 'Area';...
    'BoundingBox'; 'BoundingBox'; 'BoundingBox'; 'BoundingBox'; 'Extent';...
    'EquivalentDiameter'; 'EulerNumber'; 'Perimeter'; 'perim/area'; 'majAx/minAx';...
    '(majAx/minAx)/extent'};

%% GET folder Names
dirOutput = dir(fullfile(Alldata, 'PIC_*')); %specifiy source folder
DateFldrNms = {dirOutput.name}';
if isempty(DateFldrNms); error('No "PIC_..." folder, Check folder name and paths'); end

%% for each "PIC" folder check to see that the required files are there
% If not already done Crop the useful area from the picture
CheckGetCrop(Alldata, DateFldrNms, imgfmt)

% If filter params are missing, have user identify 5 worms and get params
Particleparams (resz, Alldata, DateFldrNms, imgfmt, dynamicTH, thresh_hold)

tic
%% loop folders
for W=1:length(DateFldrNms)
    %% settup
    RUNfinalDir = [Alldata, 'RESULTS', filesep, TrialName, filesep, DateFldrNms{W} 'RUNfinal'];
    ErrorDir = [Alldata, 'RESULTS', filesep, TrialName, filesep, DateFldrNms{W} 'ErrorDir'];
    mkdir(RUNfinalDir);
    mkdir(ErrorDir);
    %folder specicic partilce parameters
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
    
    
    DateFldrNms2 = {dirOutput2.name}'; % cell array to matrix
    CropPar=load([Alldata filesep DateFldrNms{W} filesep 'CropParam.mat']);%reload each time
    mask =CropPar.mask; posctr=CropPar.posctr;
    close all
    Images=[];
    poshead=[];
    
    for ImN=1:size(DateFldrNms2,1); % image loop
        varStruct=[] ;
        % get image and image characteristics
        startimageName=[DateFldrNms2{1}];
        imageName=[DateFldrNms2{ImN}];
        suffix=DateFldrNms{W}(findstr('_',DateFldrNms{W})+1:end);
        startImgCt=str2double(startimageName(findstr(suffix, startimageName)+length(suffix):findstr('.', startimageName)-1));
        CurrimageCt=str2double(imageName(findstr(suffix, imageName)+length(suffix):findstr('.', imageName)-1));
        imageCt=CurrimageCt-startImgCt+1;
        timeintv =imageCt/framerate;
        img=imread([Alldata filesep DateFldrNms{W} filesep DateFldrNms2{ImN}]);
        if size(img,3)>2; img=rgb2gray(img); end
        %% Remove blotchy background %works but takes a lot of time
        if strcmpi (smoothbkg, 'y')
            SE=strel('disk',10);
            imgTH=imbothat(img, SE); %works but takes a lot of time
            imgTH=imcomplement(imgTH);
            imgAdj=adapthisteq(imgTH); %Opt?
            if (strcmpi (allow_img, 'y')); figure; imshow(imgTH); figure; imshow(img1); figure; imshow(imgAdj); end
            img1=uint8(imgAdj);
        else
            img1=uint8(img);
        end
        img1=imresize(img1, resz); %resize after smoothing to save time
        close all
        
        %% MASK OFF THE PERIMITER
        mask=uint8(mask);
        img1=(img1.*mask);
        lng1=length(img1(1,:)); lng2=length(img1(:,1));
        %% Check for particles until you get a good worm or a failure
        %if masked image fails will recheck unmasked image
        particleCheck = 'in_progress'; %is switched off after second pass
        MaskImage = 'y'; %is switched off in second pass
        if exist('boundingBox')==0;  boundingBox=[1  1 size(img1,2) size(img1,1)]; end
        while strcmpi(particleCheck, 'in_progress')
            %% Area restricted search for worm uses LAST WORM BOUNDING BOX
            [img1Masked] = AdaptBoundBox(boundingBox, PadPrc, MaskImage, img1, allow_img);%allow_img
            %% MASK BY Intensity ** removing middle tones %can lead to washing out image
            close all
            if strcmpi (intenseMsk, 'y')
                [img1Masked]=IntenseMask (img1Masked, dynamicBndLim, Val, EvenImgBgSub, allow_img);
            end
            %% IDENTIFY OBJECTS and FILTER DATA
            if strcmpi(invertImage, 'y'); img1Masked=imcomplement(img1Masked); end
            if (strcmpi(dynamicTH, 'y'));
                imgBW=abs(double(img1Masked));
                thresh_hold = graythresh(imgBW);% dynamically determine threshold
                imgBW=imcomplement(im2bw(imgBW,thresh_hold));% apply threshold
            elseif (strcmpi(dynamicTH, 'y-ShortCircuit'));
                imgBW=abs(img1Masked); %SKIPS OVER INTENSITY AND DOUBLE opreations
                thresh_hold = graythresh(imgBW);% dynamically determine threshold
                imgBW=im2bw(imgBW,thresh_hold);% apply threshold
            elseif (strcmpi(dynamicTH, 'n'))
                imgBW=abs(img1Masked);%imgBW=abs(double(MasImg));
                imgBW=imcomplement(im2bw(imgBW,FltrParams.threshold));% apply threshold
            end
            
            if strcmpi(allow_img, 'y'); figure; imagesc(imgBW); end
            Image_PropertiesAll=[]; F=[]; AlabeldAll=[]; gnumAll=[]; yy=[]; yyy=[]; xx=[]; xxx=[]; ym=[];xm=[];%mask=[];
            [imgBWL, F, Image_PropertiesAll] = GetImgProps (imgBW, allow_img);
            
            %values absent - make one dummy line
            if size(Image_PropertiesAll, 1) < 1; Image_PropertiesAll= ones(1,21); end
            
            %% Filter and present data.
            Img_Propfilt=[];
            if strcmpi ('y', allow_img); figure ; imshow(imgBWL); end
            
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
            
            %% APPLY FILTERS
            countgood = 'n';
            FiltersTried=1;
            while strcmpi('n', countgood) %apply different filters to try to get a "good count"
                switch FiltApp
                    case 'all'
                        filter= [BndBxFlt4L.* BndBxFlt4U.* szFltU.* szFltL.*EccFltU.* EccFltL.*MajAxFltU.* MajAxFltL.*MinAxFltU.* MinAxFltL.*ExtFltL.*ExtFltU];
                    case 'SZ_Ax_BB'
                        filter= [szFltU.* szFltL.*MajAxFltU.* MajAxFltL.*MinAxFltU.* MinAxFltL.*TotAxFltU.*TotAxFltL];    %BndBxFlt4L.* BndBxFlt4U.*
                    case 'uin'
                        filter= filter; %specified as filter by callng program
                end
                clear ('filtVal')
                filtVal=find(filter); %finds row addresses of values
                Img_Propfilt=Image_PropertiesAll(filtVal,:); %new matrix with 0s filtered out
                
                %finally find the closest particle to the last worm identifed
                if ImN == 1; postn=FltrParams.StartPos; else postn=poshead; end
                [row,mindiff]=CloseCentr(postn,Img_Propfilt); 
                Img_Propfilt=Img_Propfilt(row, :); %select the "worm" closest to the last worm
                    % this forces a single worm thus subverting the
                    % refiltering routines below
                
                %% CLASSIFY PROBELM CASES, redo or discard image
                %too many worms, have not refiltered, then refilter.
                if and(size(Img_Propfilt, 1) > MinWorms*MaxWormFactor, FiltersTried < MaxFilt);
                    countgood = 'n';
                    StringentFilter = 'y';
                    FiltersTried=FiltersTried+1;
                    particleCheck = 'in_progress';
                 %too many worms, have refiltered, then reject.
                elseif and (size(Img_Propfilt, 1)> MinWorms*MaxWormFactor, FiltersTried == MaxFilt);
                    % if you have tried all the filters and there are still too many
                    % worms, discard the time point
                    countgood = 'n'; %skip image display and mark for continue when out of loop
                    Skiplog={Skiplog; imageName} %IS THIS POPULATED? DOES IT REALLY AVOID SAVING?
                    particleCheck = 'Done';
                    display (['done   ', num2str(ImN) 'count error']);toc
                    break % next iteration without saving in data
               elseif size(Img_Propfilt, 1)== 0;%No Particles, rerun without masking
                    switch MaskImage
                        case 'y'
                        particleCheck = 'in_progress' %is switched off after second pass
                        MaskImage = 'n';
                        countgood = 'y'; %count is not good BUT, this allow escape from REFILTERING LOOP
                        display (['unmask ', num2str(ImN)]);
                       
                        %already unmasked still no particles, REJECT
                        case'n'
                        particleCheck = '' %is switched off after second pass
                        MaskImage = 'n';
                        toc
                        display (['unmask ', num2str(ImN)]);
                        countgood = 'n'; %skip image display and mark for continue when out of loop
                        particleCheck = 'Done';
                        display (['done   ', num2str(ImN) 'ZeroCounterror']);
                        break % next iteration without saving in data
                    end
                else %OK GOOD proceede to next step
                    countgood = 'y';
                    particleCheck = 'Done';toc
                    display (['done   ', num2str(ImN)]);
                end
            end
        end
        
        %% if you have a good count then proceede
        if strcmpi('y', countgood)
            %% COLLECT DATA
            centroidLs= [centroidLs; Img_Propfilt(1,1:2)];% add a centroid
            centrSmthY=smooth(centroidLs, SmoothInt, SmoothMeth); % smooth in one set
            centrSmthY=reshape(centrSmthY, size( centroidLs, 1), size( centroidLs, 2));
            
            if length(centrSmthY) == size(centroidLs, 1)-1;
                centrSmthY=[centroidLs(1,:);centrSmthY]; % need to add a spacer.
            end
            
            CentrComp= [centroidLs, centrSmthY]; % rebuilding CentrComp each time
            Imagesfilt={};
            for ChosenIMAGE=1:length(filtVal);
                Imagesfilt = [Imagesfilt; F(filtVal(ChosenIMAGE,1)).Image];
            end
            
            %% EXTRACT THE object parameters OF INTEREST
            
            boundingBox=Img_Propfilt (1, 11:14);
            boundingBox(3)= boundingBox(3)-1;
            boundingBox(4)= boundingBox(4)-1;
            ImgCellPic=double(imcrop(img1, boundingBox));
            ImgCellmsk=double(Imagesfilt{row,1});
            ImgCell=(ImgCellmsk.*ImgCellPic);
            %% extract particle params
            areaboundingBox= boundingBox(3)* boundingBox(4);
            ImgCellNoZero=ImgCell(ImgCell~=0);
            areacell = (Img_Propfilt(1,3));
            %major/minor axis
            if boundingBox(3)> boundingBox(4) ;
                MajvsMin=boundingBox(3)/ boundingBox(4);
            else
                MajvsMin=boundingBox(4)/ boundingBox(3);
            end
            BBratio=[BBratio; MajvsMin, imageCt];
            numObj=length (Img_Propfilt (:,1));
            
            %% DISPLAY results
            if strcmpi(allow_img, 'y')
                figure; %HISTOGRAM filtered values
                title ('all filters applied')
                subplot (3,4,1); hist (Img_Propfilt  (:,5)); title ('Eccentr 1=straight')
                subplot (3,4,2); hist (Img_Propfilt  (:,8)); title ('MajAxis')
                subplot (3,4,3); hist (Img_Propfilt (:,9)); title ('MinAxis')
                subplot (3,4,4); hist (Img_Propfilt (:,10)); title ('Area')
                subplot (3,4,5); hist (Img_Propfilt  (:,15)); title ('Extent')
                subplot (3,4,6); hist (Img_Propfilt (:,16)); title ('equivalentDiameter')
                subplot (3,4,7); hist (Img_Propfilt  (:,14)); title ('EulerNumber')
                subplot (3,4,8); hist (Img_Propfilt (:,18)); title ('Perimiter/area')
                subplot (3,4,9); hist (Img_Propfilt  (:,20)); title ('maj_min_extnt')
                subplot (3,4,10); hist (Img_Propfilt  (:,14)); title ('BOundingBox4(col-14')
                %DISPLAY RESULTS
                close all
                figure; imshow(imgBW); %axis equal;
                hold on; plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'w*')
                figure ('position', scrsz); imagesc(img1); %axis equal;
                hold on; plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*')
                %ADD PARTICLE NUMBER
                for rowN = 1:length(Img_Propfilt(:,1))
                    text(Img_Propfilt(rowN,6),...
                        Img_Propfilt(rowN,7),...
                        num2str(Img_Propfilt(rowN,4)),'Color', 'r', 'FontSize', 12)
                end
                figure ('position', scrsz); imagesc(img-2); %axis equal;
                hold on; plot (Image_PropertiesAll(:,6), Image_PropertiesAll(:,7),'wO')
            end
            Image_PropertiesList = Image_PropertiesAll; ctX = 6; ctY = 7; TxtCol = 4; NumCol = 'r'; FntSz = 6;
            if strcmpi ('y', allow_img); PrintObjParamsAug11; end
            %% FILTERED FIGURE FOR SCORING DIAGNOSIS
            if strcmpi(ProofingImgs, 'y')
                nameProof= [DateFldrNms{W},'MultView',imageName,'.pdf'];
                ProofImages(scrsz, img1, Img_Propfilt, img, nameProof) %make function
            end %end proofing images loop
            %% CAPTURE AS STACKED TIFF
            if strcmpi(DataCapMode, 'StackGiff') % Append to stack of images
                Filename = [RUNfinalDir filesep DateFldrNms{W},'stack.gif'];
                measure = ['BB-ratio',num2str(MajvsMin)];%
                textls={'Thrashing in Buffer'; imageName; [num2str(timeintv), 'secs']; measure};
                Nm2=strrep(Filename, '_', '-');
                saveImageToStack(uint8(img1), Filename, ...
                    'title', 'Image', ...
                    'image_name', Nm2, ...
                    'scale_image', false, ...
                    'display_image', 'on',...
                    'CentrComp', CentrComp,...
                    'boundingBox', boundingBox,...
                    'plotcol',[1,2,3,4]); %proofingImgVIS
            end
        end
        close all
        %% COMPILE DATA, reduce image sizes
        img=imresize(img, .25); %shrink
        img1=imresize(img1, .25);
        imgBWL=imresize(imgBWL, .25);
        varStruct.filters.filtVal=filtVal;
        varStruct.images.imgBWL=imgBWL;
        varStruct.images.img1=img1; %original image
        varStruct.analysis.ImgPropHeaders=ImgPropHeaders;
        varStruct.analysis.Image_PropertiesAll=Image_PropertiesAll;
        
        %% COUNT ERROR-  Too many particles, then SKIP THE PAIR
        if strcmpi('n', countgood)
            varStruct.filters.szFltL= szFltL;
            varStruct.filters.szFltU= szFltU;
            varStruct.filters.MajAxFltL= MajAxFltL; %rows less than 1280 are in chanel 1
            varStruct.filters.MajAxFltU= MajAxFltU;% the eccentricity of the spot should be less than .8 (usually ~.5)
            varStruct.filters.MinAxFltL= MinAxFltL; %rows less than 1280 are in chanel 1
            varStruct.filters.MinAxFltU=MinAxFltU;
            varStruct.analysis.F=F;
            saveThis([ErrorDir filesep imageName(1:end-4), 'CountError.mat'], varStruct)
            continue
        end
        %% get head position for first worm
        %if there is more than 1 image choose the closest to the point
        %last head
        %>[row,mindiff]=CloseCentr(FltrParams.StartPos,Img_Propfilt);
        [WmImgPad]=GetPadImg (pad, (Imagesfilt{row,:}));
        if isempty(poshead)
            %display a few images to tell which part is the head
            Flipbook([Alldata filesep DateFldrNms{W}], DateFldrNms2(1:10));
            %function [poshead, WmImgPad]=GetHeadPosPad (pad, Imagesfilt,)
            mssg='drag point to head and double click';
            [poshead]=GetHead (poshead, WmImgPad, mssg);
            varStruct.Pos.poshead=poshead; close all;
        end
        %% spine worm
        figure; subplot(1,2,1); imshow(img1); subplot(1,2,2);  imshow(WmImgPad);
        [SpineData, poshead2] = SpineWorm (WmImgPad, allow_img, poshead, numpts);
        [poshead]=updatePoshead (poshead2, poshead);
        close all
        
        %% BADSPINE DUMP - check for errors
        if strcmpi('n', SpineData.spinegood)
            saveThis ([ErrorDir filesep imageName(1:end-4), 'SpineError.mat'],varStruct);
            continue
        end
        FltNm2={}; FltNm2= ['sze'; 'ecc'; 'mjx'; 'mnx'; 'ext'];
        FltPram2= [LowLim UpLim; eccL eccU; MajAxL MajAxU; MinAxL MinAxU; ExtL ExtU];
        
        % list of image attributes
        if strcmpi ('y', allow_img); figure; imagesc(imgBWL); end
        numObj=length (Img_Propfilt (:,1));
        Allsizes = Img_Propfilt(:,3); %look up "regionprops" in
        sizesScored = Image_PropertiesAll(:,3);
        
        %% PREPARE for separate save function (required for parallel processing)
        %Flesh out varStruct
        varStruct.images.Imagesfilt=Imagesfilt;
        varStruct.filters.FltNm2= FltNm2;
        varStruct.SpineData=SpineData;
        varStruct.analysis.Img_Propfilt=Img_Propfilt;
        SaveImNm=imageName(1:end-4);
        saveThis([RUNfinalDir filesep SaveImNm, 'final.mat'], varStruct);%'ProcessDate'
        
        %create OneWormData - for CLOSEST worm (row=1; only one row left)
        majAx_minAx = (Img_Propfilt(1,20));
        majA_minAx_extent = (Img_Propfilt(1,21));
        CentrX = (Img_Propfilt(1,1));
        CentrY = (Img_Propfilt(1,2));
        areacell = (Img_Propfilt(1,3));
        boundingBox=Img_Propfilt (1,11:14);
        OneWorm_Data = [OneWorm_Data;  {DateFldrNms{W}, ceil(numObj),...
            areacell, areaboundingBox, majAx_minAx, majA_minAx_extent,...
            CentrX, CentrY}];
        %AliveData={AbsTimeDays, timeintv, DateFldrNms{W}, numObj}
        %cellwrite(DataDir/, AliveData, 2)
    end % IMAGES LOOP
    
    OneWorm_DataHeaders =  ['filename',',', 'ObjecCount',',',...
        'areacell',',','areaboundingBox',',','majAx_minAx',',', 'majA_minAx_extent',...
        ',','CentrX',',', 'CentrY'];
    
    save ([RUNfinalDir, filesep ,DateFldrNms{W},'_',DateFldrNms2{ImN}, 'AllData.mat']) % OutDir
    save([RUNfinalDir, filesep ,DateFldrNms{W} 'OneWormdata.mat'], 'OneWorm_DataHeaders', 'OneWorm_Data')
end % FOLDERS LOOP
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


function ProofImages(scrsz, img1, Img_Propfilt, img, nameProof) %make function
%SCORED OBJECTS MULTIPLE VIEWS
figure ('position', scrsz);
subplot(2,2,1); imshow(uint8(img1))
hold on;
plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 5)
subplot(2,2,3); imagesc(img); title('subtracted images scaled '); colorbar
hold on;
plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 2)
subplot(2,2,4); imshow(img); title('subtracted images unscaled '); colorbar
hold on;
plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 2);...
    %title (['FilterNumber_',num2str(FiltersTried)])
saveas (gcf, nameProof)% save as PDF
end

function [MasImg]=IntenseMask (img, dynamicBndLim, Val, EvenImgBgSub, allow_img)
% for 'static' dynamicBndLim optio, Val is a matrix [high, low]
% for 'prc' and 'stdv' options Val is  BndLim and NumStd for
% dynamically determine bound limits from image characteristics
switch dynamicBndLim;
    case 'stdv'
        Minbnd=Val*(std2(img));  %3 std is not nearly enough
        Maxbnd=Val*(std2(img));
    case 'prc'
        Minbnd=Val*(min(min(img))); %specifically remove the middle 'bnd' % of colors
        Maxbnd=Val*(max(max(img)));
    case 'static'
        Minbnd=Val(1,2); %specifically remove the middle 'bnd' % of colors
        Maxbnd=Val(1,1);
end
avg=mean(img(:));
Mask=(img < (Minbnd) | img > (Maxbnd)); %threshold set as bounds around AVERAGE

if strcmpi(EvenImgBgSub, 'y');% only for non-simple mode
    Mask=imcomplement(Mask); % off for straight subtract ON for adjusted subtract
end
%% get STATS for the mask
imgRange=[(max(max(img))); min(min(img)); range(range(img));avg; BndLim];
RangeRmvd=[Minbnd;Maxbnd];
totPixels= size(Mask, 1) *size(Mask, 2);%Mask=(img > (avg+bnd)| img < (avg-bnd)); % perhaps substiute 0 for (avg-bnd) here?
pxlmsk=totPixels-sum(sum(Mask));
PrcMasked=(pxlmsk./totPixels).*100; %Mask=Mask-1; % masked values come in as ones so change them to zeros

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


