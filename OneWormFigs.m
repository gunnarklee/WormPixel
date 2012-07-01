function OneWormFigs(varargin)

%G. Kleemann - 10/27/11
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Version History

% Parser for case when inputdir and outputdir are specified
p = inputParser;
p.addRequired('inputdir', @isdir);
p.addRequired('outputdir', @isdir);
p.addOptional('trialname', datestr(now), @ischar);
%p.addOptional('EventList','EventList.csv', @ischar); %CSV file plate, well, AssayStart, TimeOn, TimeOff

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

%update directories to find the results file
DataDir= [DataDir 'DoneRESULTS']

%% Make new directories
ErrorDir = [Alldata, 'RESULTS', filesep, TrialName, 'ErrorDir'];
ProcessDate=(date);
SCALETYPE='(Pixels)' %unless calibrated to mm

%SummaryFiles=[DataDir RecentFldr filesep 'SummData'],

%% SPECIFY PARAMETERS
OneWorm_CHR7Params
%% get pixel/mm calibration if there is a scale image

calib= input ('Is there a claibration image? (y/n)', 's');
if strcmpi(calib, 'y');
    ScaleImage=uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'}, 'select the scale Image stretch the line to 1mM', Alldata)%InputDir
    ScaleIM=imread([Alldata filesep ScaleImage]);
    figure; imshow(ScaleIM);
    h = imline;
    position = wait(h);
    SCALETYPE='(mM)';
    %pdist2([position(1,1),position(1,2)],[position(2,1),position(2,2)])%SAME
    pixelPERmm=sqrt((position(1,1)-position(2,1))^2 + (position(1,2)-position(2,2))^2 )
end



%% load the data summary files
% Setup variables
scrsz = get(0,'ScreenSize');

%% GET THE MOST RECENT "final.mat" folder for looping
%[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'final', 'final');

centroid=[];
SpineCoords={};
CurveMtx=[];
imgs={};
CurveMtx=[];
SpineList={};
Pointlist={};
poshead=[];
velocity=[];
time=[];

%%
%update directories to find the results file

[~, RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', 'recent')

for y=1:length(namedate(:,2));%cycle folders
    
    %% GET list of "final.mat" files in folder for looping
    [DateFileNms RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', namedate(y,2));
    
    
    %% check order, sort and add spacers for skipped or missing fra
    %Namelist has spacers, DateFileNmsDoes not
    
    [NameList]= GetImgNumOrdr(DateFileNms, '-', 'final');
    
    
    %% stack the data from "final.mat" files
    %preallocate Matricies
    
    leng=length(DateFileNms);
    PaddedLeng=length(NameList);
    
    distanceMv=zeros(leng, 1);   %velocity=[velocity;vel];
    time=zeros(leng, 1);   %time=[time;insttime+time(end,:)];
    
    centroid=zeros(leng, 2);   %centroid=[centroid;CurrCent];
    imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
    CurveMtx=zeros(numpts-2, PaddedLeng);%CurveMtx=[CurveMtx,SpineData.AngleLs];
    SpineList=cell(leng, 1);   %SpineList=[SpineList,SpineData.SpineList];
    Pointlist=cell(leng, 2);   %Pointlist=[Pointlist,SpineData.Pointlist];
    
    DateFileNms_Count=0%set up separate cont for date file names
    
    for PaddedList_Count=1:length(NameList);
        PaddedList_Count
        NameList{PaddedList_Count,1}
        if NameList{PaddedList_Count,1} == 'X';
            %insert spacers in data matricies
            % distanceMv(w,1)=X;
            % centroid(w,1:2)=[X,X];   %centroid=[centroid;CurrCent];
            %imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
            CurveMtx=zeros(numpts-2, PaddedList_Count);%CurveMtx=[CurveMtx,SpineData.AngleLs];
            SpineList=cell(leng, 1);   %SpineList=[SpineList,SpineData.SpineList];
            Pointlist=cell(leng, 2);   %Pointlist=[Pointlist,SpineData.Pointlist];
            % time=zeros(leng, 1);   %time=[time;insttime+time(end,:)];
            continue;
        else
            DateFileNms_Count=DateFileNms_Count+1%set up separate cont for date file names
            load([DataDir filesep RecentFldr filesep NameList{PaddedList_Count,1}]);
        end
        
        
        img1=varStruct.images.img1;
        [CurrCent]=FindCentr(varStruct.analysis.Img_Propfilt, 'CtrMass'); %CtrMass
        insttime=1/framerate; %seconds
        
        %% use last centroid and current centroid to calcualte velocity
        if DateFileNms_Count == 1
            distanceMv(DateFileNms_Count,1)=0;
        else
            maxrow=max(find(centroid(:,1) >1));
            if isempty(maxrow); maxrow=[1]; end
            distMv=(pdist2(centroid(maxrow,:),CurrCent));%*MicroM_Pixel)/insttime;
            %mm/pixel calibration
            if exist('Calib_Scale'); distMv=distMv*Calib_Scale; end % translate distance into CM, if a scale was loaded
            distanceMv(DateFileNms_Count,1)=distMv;
        end
        
        %% cocatenate the order corrected matricies
        centroid(DateFileNms_Count,1:2)=CurrCent;
        imgs{DateFileNms_Count}=varStruct.images.imgBWL;%not sure if necc
        CurveMtx(:,DateFileNms_Count)=varStruct.SpineData.AngleLs;
        SpineList{DateFileNms_Count}=varStruct.SpineData.SpineList;
        Pointlist{DateFileNms_Count}=varStruct.SpineData.Pointlist;
        
        %calculate times
        if DateFileNms_Count==1
            time(DateFileNms_Count)=0;
        else
            time(DateFileNms_Count)=insttime+time(DateFileNms_Count-1);
        end
        
        %% Visualise data assembly PROGRESS
        
        if strcmpi (allow_img, 'y')
            close all
            figure;imagesc(CurveMtx);
            figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,size(img1,1)]); ylim([1, size(img1, 2)]);
            %>>  figure;plot(time,velocity,'-r') ; title ('displacement vs. time');
            hold on
            axes('position', [.75 .15 .15 .15]);
            plot(SpineList{w,1}(1,1), SpineList{w,1}(1,2), 'r');
            
            %%CHECK CURVE MATRICIES
            CurveMtxtemp(:,1:15)=repmat(varStruct.SpineData.AngleLs, 1, 15);
            pos=[50 50 size(varStruct.images.Imagesfilt{1,1},1) size(varStruct.images.Imagesfilt{1,1},2)];
            
            figure; subplot(1,2,1); imagesc(CurveMtxtemp);
            [WmImgPad] = padImg (varStruct.images.Imagesfilt{1,1}, pad);
            subplot(1,2,2); imagesc(WmImgPad); colorbar;
            
            figure; imagesc(CurveMtxtemp);colorbar;
            
            
            %    WmImgPadcolor=(imoverlay (mat2gray(WmImgPad), varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), [1,0,0]));
            figure; imshow(WmImgPad, 'InitialMagnification', 400);
            hold on
            plot( varStruct.SpineData.SpineList(:,2), varStruct.SpineData.SpineList(:,1), 'r.', 'MarkerSize', 10)
            plot(Pointlist{1,1}(:,2), Pointlist{1,1}(:,1), 'b*', 'MarkerSize', 5)
            varStruct.SpineData.AngleLs
        end
        
        if strcmpi (stoppoint, 'y')
            stopPt= input ('next image?', 's');
        end
        
        close all
    end
end
%% SAVE EVERYTHING

outputDir=[Alldata filesep 'ResultsFiles']
mkdir(outputDir);
NameOut=[outputDir filesep RecentFldr(1:end-8)];

save ([NameOut 'SummData'], 'centroid', 'imgs', 'CurveMtx', 'SpineList', 'Pointlist', 'distanceMv', 'time', 'img1'); %does not save sub structure
TrialData=dataset(centroid (:,1), centroid (:,2), distanceMv, time, 'VarNames', {'X','Y','Dist','time (sec)'});
%export(TrialData, 'XLSfile', [DataDir RecentFldr 'movedata']);
%export(TrialData, 'XLSfile', [NameOut '_movedata']);
export(TrialData,'file', [NameOut '_movedata'],'Delimiter',',')

CurveData=dataset(time, CurveMtx');
%export(CurveData, 'XLSfile', [DataDir RecentFldr 'curvedata']);
%export(TrialData, 'XLSfile', [NameOut 'curvedata']);
export(TrialData,'file', [NameOut 'curvedata'],'Delimiter',',')


%% PLOT FIGURES
close all

numFr =50; % frame count for the zoomed in view
numFr2 =25;
clims = [-35 35];

%% curve matricies
%full color full scale
figure;imagesc(CurveMtx); colorbar
saveas (gcf, [NameOut 'CurveMtxrfull'], 'pdf')

MakeTricolMap

% tricolor zoomed
figure;imagesc(CurveMtx(:,1:numFr));
colormap(RWBcMap2); colorbar;
saveas (gcf, [NameOut 'CurveMtx' num2str(numFr)], 'pdf');

% tricolor all frames
figure;imagesc((CurveMtx), clims);
colormap(RWBcMap2); colorbar
%>> saveas (gcf, [NameOut 'CurveMtxallframes_limited'], 'pdf')

% tricolor no limits
figure;imagesc(CurveMtx);
colormap(RWBcMap2); colorbar
title ('Tricolor Curve Matrix allframes'); xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
saveas (gcf, [NameOut 'CurveMtxallframes'], 'pdf')

%% Path traveled
figure;plot(centroid(3:end,1),centroid(3:end,2), '-b',  'MarkerSize', 2);...
    title ('centroid position'); xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
saveas (gcf, [NameOut 'PathTraveled'], 'pdf')

%% plot displacement
figure;plot(time,distanceMv,'-r') ; title ('displacement between frames');
hold on
xlabel ('time (Sec)');
ylabel (['displacement ', SCALETYPE]);
%>>axes('position', [.75 .15 .15 .15]);
%>>plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
saveas (gcf, [NameOut 'Distvstm'], 'pdf')

%% Cummulative Dist
CumulDist=cumsum(distanceMv)
figure;plot(time,CumulDist,'-r') ; title (['Cumlative displacement' SCALETYPE ' vs. time']);
xlabel ('time (Sec)');
ylabel (['Cumulative displacement ', SCALETYPE]);

%>> hold on
%>> axes('position', [.75 .15 .15 .15]);
%>> plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
saveas (gcf, [NameOut 'Cuml_Travel'], 'pdf')

%% plot neck movement


ZeroLine=[repmat(0,numFr,1),time(1:numFr)]
OscilAngle=CurveMtx(3,1:numFr)

%CYCLE COUNT by looking for crossing of the midline (angle changes from
% + to -
[ChangeRow] = SignChangeDetect(OscilAngle);
cyclessStart = start(1:2:length(ChangeRow))
%cyclesEndsROWS=start(2:end)



figure;
plot(time(1:numFr),OscilAngle,'-b');
hold on
plot(ZeroLine(:,2),ZeroLine(:,1),'-k')
plot(cyclesEndsROWS, 0, '*r')

title ('Head oscillations')
xlabel ('time (Sec)');
ylabel (['degrees from the midline']);

saveas (gcf, [NameOut 'Neck Movement' num2str(numFr)], 'pdf')
%% PLOT CYCLE (ENDS)



%% CYCLES PER SECOND (Varn for crOss strain comparisons)



%% distance , cycle (efffort/ work ratio)


figure; plot(CurveMtx(3,:));
title ('Head oscillations')
xlabel ('time (Sec)');
ylabel (['degrees from the midline']);

saveas (gcf, [NameOut 'All Neck Movement'], 'pdf')


%     plot(SpineList{W}(:,1), SpineList{W}(:,2));
%     hold on
%     %figure; imagesc(CurveMtx);
%

%[wkSpnX,wkspnY]=ind2sub(size(SpineData.SpineList), find(SpineData.SpineList));
%cmap=colormap(jet (size(wkSpnX,1))) %copper
%WmImgPadcolor=mat2gray(WmImgPad)

%WmImgPadcolor=(imoverlay (mat2gray(img1), SpineData.SpineList,  cmap(Pt,:)));

%     figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);
%     hold on
%     %addpoint to imgage with new color
%     WmImgPadcolor=(imoverlay (mat2gray(WmImgPadcolor), skeleEND, cmap(Pt,:)));
end

function [CurrCent]=FindCentr(Img_Propfilt, mode) %CtrMass

switch mode
    case 'CtrMass'
        %centroid as center of mass of particle
        CurrCent=[Img_Propfilt(1,1),Img_Propfilt(1,2)];
        
    case 'CtrBB'
        %centroid as center of bounding box (Pierce-Shimomoura method)
        uplfX=Img_Propfilt(11);
        uplfY=Img_Propfilt(12);
        Xwid=Img_Propfilt(13);
        yhgt=Img_Propfilt(14);
        CurrCent=[(uplfX+(Xwid/2)),(uplfY-(yhgt/2))];
end
end

function [NameLst]= GetImgNumOrdr(DateFileNms, Prfx, Sufx)
%extract image number
%get proper image order and iser spacers for missing  images

NumLs=[] %get mage numbers from name
for im=1:length(DateFileNms)
    ImgNum=str2num(DateFileNms{im,1}(max(strfind(DateFileNms{im,1}, Prfx)+1):strfind(DateFileNms{im,1}, Sufx)-1));
    NumLs=[NumLs; ImgNum];
end

%generate full list
%%FullLS=(1:max(NumLs))';

%for each item on list find a row match
%[B,IX] = sort(NumLs);
%make a new list Plaing either 'X' or a image name if it exists

NameLst={};  %make a new list with spacers
for LsItm=(1:max(NumLs));
    NameLst=[NameLst;'X'];
end

for Nm=1:size(DateFileNms,1)
    NameLst(NumLs(Nm))= DateFileNms(Nm);
end
end

function [ChangeRow2] = SignChangeDetect(angleList)
%Input a list of angles and record the rows where the sign of an angle
%changes from + to - or to 0
ChangeRow=[];
signlist=sign(angleList);
for CurrSign=1:length(signlist)-1;
    if signlist(CurrSign)==signlist(CurrSign+1);
        change=0
    else
        change=1
    end
    ChangeRow=[ChangeRow;change];
end
[ChangeRow, signlist(1:end-1)'] %check
ChangeRow1=find(ChangeRow)
%filter out doubles. These are usually transiting through zero
ChangeRow2=[];
for CheckDbl=1:length(ChangeRow)-1
    if ChangeRow(CheckDbl)+ChangeRow(CheckDbl+1) == 2
    continue %skip the first of double ones
    elseif ChangeRow(CheckDbl)==1;
        ChangeRow2=[ChangeRow2;CheckDbl];%capture the row number
    end
end
end

