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

%% Make new directories
ErrorDir = [Alldata, 'RESULTS', filesep, TrialName, 'ErrorDir'];
ProcessDate=(date);
%SummaryFiles=[DataDir RecentFldr filesep 'SummData'],

%% SPECIFY PARAMETERS
OneWorm_CHR7Params

%% load the data summary files
% Setup variables
scrsz = get(0,'ScreenSize');

%% Make new directories
[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'final', 'final')

centroid=[];
SpineCoords={};
CurveMtx=[];
imgs=[];
CurveMtx=[]
SpineList={}
Pointlist={}
poshead=[];
velocity=[];
time=[];

%% check each folder or crop params if absent get them

for w=1:length(DateFileNms)
    load([DataDir 'RESULTS' filesep RecentFldr filesep DateFileNms{w,1}])
    CurrCent=[Img_Propfilt(1,1),Img_Propfilt(1,2)]
    insttime=1/framerate %seconds
    
    %use last centroid and current centroid to calcualte velocity
    if isempty(velocity)
        velocity=[0];
    else
        vel=(pdist2(centroid(end,:),CurrCent)*MicroM_Pixel)/insttime;
        %sqrt((CurrCent(2)-centroid(end,2))^2+(CurrCent(1)-centroid(end,1))^2)
        velocity=[velocity;vel]
    end
    
    
    
    %cocatenate the order corrected matricies
    centroid=[centroid;CurrCent];
    imgs=[imgs;imgBWL];
    CurveMtx=[CurveMtx,SpineData.AngleLs];
    SpineList=[SpineList,SpineData.SpineList];
    Pointlist=[Pointlist,SpineData.Pointlist];
   
   
    
    
    if isempty(time)
        time=[0]
    else
        time=[time;insttime+time(end,:)]
    end
    
    
    
%% make figs    
    figure;imagesc(CurveMtx)
    figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,size(img1,1)]); ylim([1, size(img1, 2)]); 
    figure;plot(time,velocity,'-r') ; title ('velocity vs. time')
    hold on
    axes('position', [.75 .15 .15 .15]) 
    plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r')
    
    if strcmpi (stoppoint, 'y')
        stopPt= input ('next image?', 's')
    end
    
    close all
end
  save ([DataDir RecentFldr 'SummData'], 'centroid', 'imgs', 'CurveMtx', 'SpineList', 'Pointlist', 'velocity', 'time'); %does not save sub structure
    
    figure;imagesc(CurveMtx)
    figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,size(img1,1)]); ylim([1, size(img1, 2)]); 
    figure;plot(time,velocity,'-r') ; title ('velocity vs. time')
    hold on
    axes('position', [.75 .15 .15 .15]) 
    plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r')
    
    
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
