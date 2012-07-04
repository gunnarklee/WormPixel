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
pixelPERmm=1 %default
%SummaryFiles=[DataDir RecentFldr filesep 'SummData'],

%% SPECIFY PARAMETERS
OneWorm_CHR7Params
%% get OIXEL/MM CALIBRATION
calib= input ('Is there a claibration image? (y/n)', 's');
if strcmpi(calib, 'y');
    ScaleImage=uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'}, 'select the scale Image stretch the line to 1mM', Alldata)%InputDir
    ScaleIM=imread([Alldata filesep ScaleImage]);
    ScaleIM=imresize(ScaleIM, resz);
    figure; imshow(ScaleIM);
    h = imline;
    position = wait(h);
    SCALETYPE='(mM)';
    %pdist2([position(1,1),position(1,2)],[position(2,1),position(2,2)])%SAME
    pixelPERmm=sqrt((position(1,1)-position(2,1))^2 + (position(1,2)-position(2,2))^2 )
end

%% GET THE MOST RECENT "final.mat" folder for looping
%[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'final', 'final');
scrsz = get(0,'ScreenSize');

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

%% update directories to find the results file
[~, RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', 'recent')



for y=1:length(namedate(:,2));%cycle folders
    %% GET list of "final.mat" files in folder for looping
    [DateFileNms RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', namedate(y,2));
    
    %% check order, sort and add spacers for skipped or missing fra
    %Namelist has spacers, DateFileNmsDoes not
    [NameList]= GetImgNumOrdr(DateFileNms, '-', 'final');
    
    %% stack the data from "final.mat" files
    
    insttime=1/framerate; %seconds
    
    %preallocate Matricies
    leng=length(DateFileNms);  %UnPadded count
    PaddedLeng=length(NameList); %Padded count
    distanceMv=zeros(PaddedLeng, 1);   %velocity=[velocity;vel];
    time=zeros(PaddedLeng, 1);   %time=[time;insttime+time(end,:)];
    centroid=zeros(PaddedLeng, 2);   %centroid=[centroid;CurrCent];
    imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
    CurveMtx=zeros(numpts-2, PaddedLeng);%CurveMtx=[CurveMtx,SpineData.AngleLs];
    SpineList=cell(PaddedLeng, 1);   %SpineList=[SpineList,SpineData.SpineList];
    Pointlist=cell(PaddedLeng, 2);   %Pointlist=[Pointlist,SpineData.Pointlist];
    DateFileNms_Count=0%set up separate cont for UNPADEDDED LIST date file names
    
    for PaddedList_Count=1:length(NameList); %Put the data together
        PaddedList_Count
        NameList{PaddedList_Count,1}
        
        %for things that refer to the number above
        if PaddedList_Count == 1
            distanceMv(PaddedList_Count,1)=0;
            %CurveMtx(:,PaddedList_Count)=zeros(numpts-2, PaddedList_Count);
            %SpineList{PaddedList_Count}=zeros(cell(leng, 1));
            %Pointlist{PaddedList_Count}=zeros(cell(leng, 2));
            time(PaddedList_Count)=0; %determined by FPS rate
            %centroid(PaddedList_Count,1:2)=[centroid(PaddedList_Count-1,1:2)]; %just take the last value
            DateFileNms_Count=DateFileNms_Count+1%set up separate cont for date file names
            continue
        end
        
        if and(PaddedList_Count > 1, NameList{PaddedList_Count,1} == 'X'); %then fill a blank in some matricies
            %insert spacers in data matricies
            distanceMv(PaddedList_Count,1)=distanceMv(PaddedList_Count-1,1);
            %CurveMtx(:,PaddedList_Count)=zeros(numpts-2); % JUST LEAVE IT
            %AS ZEROS
            %SpineList{PaddedList_Count}=zeros(cell(leng, 1));
            %Pointlist{PaddedList_Count}=zeros(cell(leng, 2));
            time(PaddedList_Count) = insttime+time(end,:) %determined by FPS rate
            centroid(PaddedList_Count,1:2)=[centroid(PaddedList_Count-1,1:2)]; %just take the last value
            %imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
            continue;
        end
        
        
        DateFileNms_Count=DateFileNms_Count+1%set up separate cont for date file names
        load([DataDir filesep RecentFldr filesep NameList{PaddedList_Count,1}]);
        
        
        img1=varStruct.images.img1;
        [CurrCent]=FindCentr(varStruct.analysis.Img_Propfilt, 'CtrMass'); %CtrMass
        
        %% cocatenate the order corrected matricies
        centroid(DateFileNms_Count,1:2)=CurrCent;
        imgs{DateFileNms_Count}=varStruct.images.imgBWL;%not sure if necc
        CurveMtx(:,DateFileNms_Count)=varStruct.SpineData.AngleLs;
        SpineList{DateFileNms_Count}=varStruct.SpineData.SpineList;
        Pointlist{DateFileNms_Count}=varStruct.SpineData.Pointlist;
        
        %% DISTANCE AND VELOCITY
        %use last centroid and current centroid to calcualte velocity
        % the odd loop is used because of discreopancies btwn padded and
        % unpadded list
        if DateFileNms_Count == 1
            distanceMv(DateFileNms_Count,1)=0;
        else
            LastCentr=centroid((max(find(centroid(:,1) >1)))-1,:) %find the last centroid with a value
            if LastCentr == 0, distMv = 0
            else
                distMv=(pdist2(LastCentr,CurrCent))
            end
            distanceMv(DateFileNms_Count,1)=distMv
        end
        
        time(DateFileNms_Count)=insttime+time(DateFileNms_Count-1); %calculate times
        
        %%PLACE VUSIALIZE FUNCITON HERE
        
        close all
    end
end

%% GET DERIVED PARAMETERS
%these include path traveled (smoothed displacemnt)
%and reset time for head subcycles (passed 0 degrees)

% PATH PARAMETERS
%Get the non zero values for centroid and corresonding time
centroidPlot=centroid((find(~centroid(:,1)==0)), 1:2);

% USE THIS FOR COLOR EVOLUTION
timePlot=time(find(~centroid(:,1)==0));
%[wkSpnX,wkspnY]=ind2sub(size(SpineData.SpineList), find(SpineData.SpineList));
%cmap=colormap(jet (size(wkSpnX,1))) %copper

%% smooth path with sliding window
%STILL NEEDS WORK (run of ZEROS)
smoothWin=25;
StraightWin=50;
%Smooth forward
pathSmooth=[SmoothData(centroidPlot(:,1), smoothWin), SmoothData(centroidPlot(:,2), smoothWin)];
pathStrgt=[SmoothData(centroidPlot(:,1), StraightWin), SmoothData(centroidPlot(:,2), StraightWin)];

% % take last data point and draw a line to it with intermediate points
% between
% find the point where the movemb
Breakpt=length(pathStrgt)-StraightWin*.1
pathStrgtStart=pathStrgt(1:Breakpt,1:2)
%FILL IN  W/unsmoothed data rather than the single mean position (no
%movement at all)
pathStrgtEND=[]
Startval=centroidPlot(Breakpt+1,1:2)
Endval=centroidPlot(end,1:2)
diffVal=Endval-Startval
intervals= length(centroidPlot)-(Breakpt+1)
subdiv=diffVal/intervals;
pathStrgtEND=[[Startval(1):subdiv(1):Endval(1)]', [Startval(2):subdiv(2):Endval(2)]']

%Stitch together
pathStrgt=[pathStrgtStart;pathStrgtEND]

%%
figure;plot(centroidPlot(2:end,1),centroidPlot(2:end,2), '-b',  'MarkerSize', 2);...
    hold on
%plot(pathSmooth(:,1),pathSmooth(:,2), '--r',  'MarkerSize', 2);...
  plot(pathStrgt(:,1),pathStrgt(:,2), '--k',  'MarkerSize', 2, 'lineWidth', 3);...
%
%    plot(pathStrgtStart(:,1),pathStrgtStart(:,2), '--b',  'MarkerSize', 2, 'lineWidth', 3);...
   % plot(pathStrgtEND(:,1),pathStrgtEND(:,2), '--o',  'MarkerSize', 2, 'lineWidth', 3);...
pathSmooth(length(pathStrgt)-StraightWin+1:end,1)=centroidPlot(length(pathStrgt)-StraightWin+1:end, 1)
%pathSmooth(length(pathStrgt)-StraightWin+1:end,2)=centroidPlot(length(pathStrgt)-StraightWin+1:end, 2)


%% TRAVEL Distance from STRAIGHTened path.

TravelDist=zeros(size(pathStrgt, 1),1);   %velocity=[velocity;vel];
for WmPos=2:size(pathStrgt, 1); %Put the data together
    WmPos
    distmoved=pdist2(pathStrgt(WmPos-1,:), pathStrgt(WmPos,:))
    TravelDist(WmPos)=distmoved;
end

%%

%adjust to calibrated mm/pixel (1 if not calibrated)
TravelDist = TravelDist/pixelPERmm
distanceMv = distanceMv/pixelPERmm

%% SWIM STROKE statistics
ZeroLine=[repmat(0,length(time),1),time]
OscilAngle=CurveMtx(3,1:end)'

%STROKE COUNT
%looking for crossing of the midline (angle changes from + to -
%** has to be sustained on one side for more than 1 timepoint**
[StrokeMk,  signlist] = SignChangeDetect(OscilAngle, time, insttime);


StrokeCt=length(StrokeMk)-1
DistPerStroke=sum(TravelDist)/StrokeCt

%cyclessStart = start(1:2:length(ChangeRow))
%cyclesEndsROWS=start(2:end)

%preallocate matricies
InterStrokeTime=zeros(length(StrokeMk),1)
AmplStroke=zeros(length(StrokeMk),1)

for CurrStk=2:length(StrokeMk);
    
    StrokeStartTime=StrokeMk(CurrStk-1);
    StrokeEndTime=StrokeMk(CurrStk);
    StrokeStartROW=min(find(time >= StrokeMk(CurrStk-1)));
    StrokeEndROW=min(find(time >= StrokeMk(CurrStk)));
    
    InterStrokeTime(CurrStk)=StrokeEndTime-StrokeStartTime
    
    %Max angle calc
    if signlist(CurrStk) > 0
        %StrokePolarity(CurrStk)=1;
        AmplStroke(CurrStk)=max(OscilAngle(StrokeStartROW:StrokeEndROW));
    elseif signlist(CurrStk) < 0
        AmplStroke(CurrStk)=min(OscilAngle(StrokeStartROW:StrokeEndROW));
        %StrokePolarity(CurrStk)=-1
        continue
    end
end

%% SAVE EVERYTHING
outputDir=[Alldata filesep 'ResultsFiles']
mkdir(outputDir);
NameOut=[outputDir filesep RecentFldr(1:end-8) filesep];
mkdir(NameOut);

%MOVEMENT DATA
save ([NameOut 'SummData'], 'centroid', 'imgs', 'CurveMtx',...
    'SpineList', 'Pointlist', 'distanceMv', 'time', 'img1',...
    'InterStrokeTime', 'DistPerStroke',...
    'AmplStroke', 'StrokeMk', 'signlist') %does not save sub structure
%TrialData=dataset(centroid (:,1), centroid (:,2), distanceMv, time,..
%'VarNames', {'X','Y','Dist','time (sec)'});
TrialData=dataset(centroid (:,1), centroid (:,2), distanceMv, time,...
    'VarNames', {'X','Y','Dist','time (sec)'});
%export(TrialData, 'XLSfile', [DataDir RecentFldr 'movedata']);
if ispc
    export(TrialData, 'XLSfile', [NameOut 'movedata']);
else
    export(TrialData,'file', [NameOut 'movedata'],'Delimiter',',');
end

%STROKE DATA
Strokedata=dataset(InterStrokeTime, AmplStroke, StrokeMk,...
    signlist, 'VarNames', {'StrokeDurn', 'StrokeAmpl',...
    'StrokeMk', 'signlist'}); %maybe add stroke specific power...'Strokeower'

%export(TrialData, 'XLSfile', [DataDir RecentFldr 'movedata']);
if ispc
    export(Strokedata, 'XLSfile', [NameOut 'SwimStrokedata']);
else
    export(Strokedata,'file', [NameOut 'SwimStrokedata'],'Delimiter',',');
end

%CURVE DATA
CurveData=dataset(time, CurveMtx');
%export(CurveData, 'XLSfile', [DataDir RecentFldr 'curvedata']);
if ispc
    export(CurveData, 'XLSfile', [NameOut 'curvedata']);
else
    export(CurveData,'file', [NameOut 'curvedata'],'Delimiter',',')
end

%% PLOT FIGURES
close all

numFr =50; % frame count for the zoomed in view
%>numFr2 =25;
%>clims = [-35 35];

%% curve matricies
%full color full scale
figure;imagesc(CurveMtx); colorbar
%saveas (gcf, [NameOut 'CurveMtxrfull'], 'pdf')

MakeTricolMap

% tricolor zoomed no limits
figure;imagesc(CurveMtx(:,1:numFr));
colormap(RWBcMap2); colorbar;
saveas (gcf, [NameOut 'CurveMtx' num2str(numFr)], 'pdf');

% tricolor all frames
%>figure;imagesc((CurveMtx), clims);
%>colormap(RWBcMap2); colorbar
%>> saveas (gcf, [NameOut 'CurveMtxallframes_limited'], 'pdf')

% tricolor no limits all frames
figure;imagesc(CurveMtx);
colormap(RWBcMap2); colorbar
title ('Tricolor Curve Matrix allframes');% xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
saveas (gcf, [NameOut 'CurveMtxallframes'], 'pdf')

%% Path traveled
figure;plot(centroidPlot(2:end,1),centroidPlot(2:end,2), '-b',  'MarkerSize', 2);...
    title ('centroid position'); 
    xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
saveas (gcf, [NameOut 'centroid location per frame'], 'pdf')

%% PLOT Smooth path -
figure;plot(centroidPlot(2:end,1),centroidPlot(2:end,2), '-b',  'MarkerSize', 2);...
    hold on
    %plot(pathSmooth(:,1),pathSmooth(:,2), '--r',  'MarkerSize', 2);...
    plot(pathStrgt(:,1),pathStrgt(:,2), '--k',  'MarkerSize', 2, 'lineWidth', 3);...
    legend ('centroid path', 'path traveled (smoothed)')
    title ('travel path');
%xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
saveas (gcf, [NameOut 'PathTraveled'], 'pdf')


%% plot displacement

figure;plot(time(1:end-1),distanceMv(1:end-1),'-r') ; title ('Centroid displacement between frames');
%last frame is zero,
hold on
xlabel ('time (Sec)');
ylabel (['centorid movement', SCALETYPE]);
%>>axes('position', [.75 .15 .15 .15]);
%>>plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
%>>saveas (gcf, [NameOut 'centroid movement'], 'pdf')


%% plot PATH TRAVELED  displacement
figure;plot(time(2:end-1), TravelDist(1:end),'-r') ; title ('displacement between frames');
%last frame is zero,
%hold on
xlabel ('time (Sec)');
ylabel (['centorid movement', SCALETYPE]);
%>>axes('position', [.75 .15 .15 .15]);
%>>plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
%saveas (gcf, [NameOut 'DIstance Traveled'], 'pdf')

%% Cummulative Dist
CumulDist=cumsum(TravelDist);
figure;plot(time(1:length(CumulDist)),CumulDist,'-r');
title (['Cumlative travel (straight path)' SCALETYPE ' vs. time']);
xlabel ('time (Sec)');
ylabel (['Cumulative Travel distance ', SCALETYPE]);

%>> hold on
%>> axes('position', [.75 .15 .15 .15]);
%>> plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
saveas (gcf, [NameOut 'Cuml_Travel'], 'pdf')

%% plot SUBSET neck movement 
figure;
plot(time(1:numFr),OscilAngle(1:numFr),'-b');
hold on
plot(ZeroLine((1:numFr),2),ZeroLine((1:numFr),1),'-k')
%plot(StrokeMk, 0, '*r')

title (['Head oscillations' num2str(numFr) 'frames'])
xlabel ('time (Sec)');
ylabel ('degrees from the midline');

saveas (gcf, [NameOut 'Neck Movement' num2str(numFr)], 'pdf');
%% plot FULL RUN neck movement 
fullRun=length(time)-1 %puts zeros at tend 

figure;
plot(time(1:fullRun),OscilAngle(1:fullRun),'-b');
hold on
plot(ZeroLine((1:fullRun),2),ZeroLine((1:fullRun),1),'-k')
plot(StrokeMk, 0, '*r'); 
h=legend ('head angle','stroke reset line', 'stroke marker')
set (h, 'Location', 'SouthEast');


title ('Head oscillations Full Run')
xlabel ('time (Sec)');
ylabel (['degrees from the midline']);
xlim([0,max(time)]);

saveas (gcf, [NameOut 'Neck MovementFullRun'], 'pdf')
%% PLOT CYCLE (ENDS)



%% CYCLES PER SECOND (Varn for crOss strain comparisons)



%% distance , cycle (efffort/ work ratio)


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

close all
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

NumLs=[]; %get mage numbers from name
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

function [StrokeZero, signlist2] = SignChangeDetect(angleList, time, insttime)
%Input a list of angles and record the rows where the sign of an angle
%changes from + to - or to 0
ChangeRow=[];
ChangeRow2=[];
signlist2=[];
ZeroCheck=[]; %Is the angle zero?
StrokeZero=[];

signlist=sign(angleList); % + or -
for CurrSign=1:length(signlist)-1; % mark location of sign changes
    if signlist(CurrSign)==signlist(CurrSign+1); %if sign of angle stays same
        change=0;
    else %if sign of angle changes
        change=1;
    end
    ChangeRow=[ChangeRow;change];
end
%[ChangeRow, signlist(1:end-1)'] %check
ChangeRow1=find(ChangeRow); % get list of rows with changes


%% FILTER OUT DOUBLES CONCATENATE DATA
%filter out doubles.These are usually transiting through zero
for CheckDbl=1:length(ChangeRow)-1
    if ChangeRow(CheckDbl)+ChangeRow(CheckDbl+1) == 2 %skip double
        continue
    elseif ChangeRow(CheckDbl)==1; % single cahnge capture data
        if angleList(CheckDbl) == 0
            StrokeZero=[StrokeZero;time(CheckDbl)] % capture the time
        else % head does not cross at measured time capture time + 1/2 iterval
            StrokeZero=[StrokeZero;time(CheckDbl)+insttime/2] % capture the time
        end
        
        %capture the sign
        if CheckDbl == 1
            signlist2=[signlist2;0];
        else
            signlist2=[signlist2;signlist(CheckDbl-1)]; %capture the sign immediately before head crosses zero
        end
    end
end
%StrokeZero=time(ChangeRow2)
end

function visualize
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

end

function countSm=SmoothData(Count, winSz)% smooth with sliding window
%% 10/2011 G. kleemann
countSm=Count;%set up array to be smoothed
for ColNm=1:size(Count,2);
    for RowNm=1:length(Count(:,ColNm))-winSz;
        TargetData=Count((RowNm:RowNm+winSz-1),ColNm);%display Target data in window
        medCt=median(Count((RowNm:RowNm+winSz-1),ColNm)) % get median for window
        replPosn=RowNm+ceil(winSz/2)
        RowNm
        countSm(replPosn ,ColNm)=medCt; % assign median to
    end
    %fill end with median
    countSm(length(Count(:,ColNm))-winSz+1:length(Count(:,ColNm)))=medCt
end
end