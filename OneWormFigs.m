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

%% GET THE MOST RECENT "final.mat" folder for looping
[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'final', 'final');

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


%preallocate Matricies

leng=length(DateFileNms)

distanceMv=zeros(leng, 1);   %velocity=[velocity;vel];
centroid=zeros(leng, 2);   %centroid=[centroid;CurrCent];
imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
CurveMtx=zeros(numpts-2, leng);%CurveMtx=[CurveMtx,SpineData.AngleLs];
SpineList=cell(leng, 1);   %SpineList=[SpineList,SpineData.SpineList];
Pointlist=cell(leng, 2);   %Pointlist=[Pointlist,SpineData.Pointlist];
time=zeros(leng, 1);   %time=[time;insttime+time(end,:)];


%% stack the data from "final.mat" files

for w=1:length(DateFileNms);
    load([DataDir filesep 'RESULTS' filesep RecentFldr filesep DateFileNms{w,1}]);
    
    [CurrCent]=FindCentr(Img_Propfilt, 'CtrMass'); %CtrMass
    
    insttime=1/framerate; %seconds
    
    %use last centroid and current centroid to calcualte velocity
    if w==1
        distanceMv(w,1)=0;
    else
        distMv=(pdist2(centroid(end-1,:),CurrCent));%*MicroM_Pixel)/insttime;
        distanceMv(w,1)=distMv;
    end
    
    %cocatenate the order corrected matricies
    centroid(w,1:2)=CurrCent;
    imgs{w}=imgBWL;%not sure if necc
    CurveMtx(:,w)=SpineData.AngleLs;
    SpineList{w}=SpineData.SpineList;
    Pointlist{w}=SpineData.Pointlist;
    
    if w==1
        time(w)=0;
    else
        time(w)=insttime+time(w-1);
    end
    
    
    
    %% make figs
    
    if strcmpi (allow_img, 'y')
        figure;imagesc(CurveMtx);
        figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,size(img1,1)]); ylim([1, size(img1, 2)]);
      %>>  figure;plot(time,velocity,'-r') ; title ('displacement vs. time');
        hold on
        axes('position', [.75 .15 .15 .15]);
        plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r');
        
        %%CHECK CURVE MATRICIES
        CurveMtxtemp(:,1:15)=repmat(SpineData.AngleLs, 1, 15);
        pos=[50 50 size(Imagesfilt{1,1},1) size(Imagesfilt{1,1},2)]

        figure; subplot(1,2,1); imagesc(CurveMtxtemp)
        [WmImgPad] = padImg (Imagesfilt{1,1}, pad)
        subplot(1,2,2); imagesc(WmImgPad)
        
        WmImgPadcolor=(imoverlay (mat2gray(WmImgPadcolor), skeleEND,  cmap(Pt,:)));
        figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);
        hold on
        plot(Pointlist(:,2), Pointlist(:,1), 'b*', 'MarkerSize', 10)
        
        
    end
    
    if strcmpi (stoppoint, 'y')
        stopPt= input ('next image?', 's');
    end
    
    close all
end


numFr =250

save ([DataDir RecentFldr 'SummData'], 'centroid', 'imgs', 'CurveMtx', 'SpineList', 'Pointlist', 'distanceMv', 'time', 'img1'); %does not save sub structure

MakeTricolMap
figure;imagesc(CurveMtx)
colormap(RWBcMap2); colorbar

saveas (gcf, [DataDir RecentFldr 'CurveMtx'], 'pdf')

figure;imagesc(CurveMtx(:,1:numFr))
colormap(RWBcMap2); colorbar
saveas (gcf, [DataDir RecentFldr 'CurveMtx' num2str(numFr)], 'pdf')


figure;plot(centroid(:,1),centroid(:,2), '--rs',...
    'LineWidth', 2, 'MarkerEdgeColor','k',...
    'MarkerFaceColor','g', 'MarkerSize', 1);...
    title ('centroid position'); xlim ([1,size(img1,2)]); ylim([1, size(img1, 1)]);
saveas (gcf, [DataDir RecentFldr 'PathTraveled'], 'pdf')

figure;imagesc(CurveMtx)
figure;plot(centroid(1:numFr,1),centroid(1:numFr,2), '--rs',...
    'LineWidth', 2, 'MarkerEdgeColor','k',...
    'MarkerFaceColor','g', 'MarkerSize', 1);...
    title ('centroid position'); xlim ([1,size(img1,2)]); ylim([1, size(img1, 1)]);
saveas (gcf, [DataDir RecentFldr 'PathTraveled' num2str(numFr)], 'pdf')



%figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,450]); ylim([1,400]);

%plot displacement
figure;plot(time,distanceMv,'-r') ; title ('displacement vs. time');
hold on
axes('position', [.75 .15 .15 .15]);
plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r');
saveas (gcf, [DataDir RecentFldr 'Distvstm'], 'pdf')

figure;plot(time(1:numFr,:),distanceMv(1:numFr,:),'-r') ; title ('displacement vs. time');
hold on
axes('position', [.75 .15 .15 .15]);
plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r');
saveas (gcf, [DataDir RecentFldr 'DistVstm' num2str(numFr)], 'pdf')

%Cummulative Dist
CumulDist=cumsum(distanceMv)
figure;plot(time,CumulDist,'-r') ; title ('Cumldispl vs. time');
hold on
axes('position', [.75 .15 .15 .15]);
plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r');
saveas (gcf, [DataDir RecentFldr 'CumlPathTraveled'], 'pdf')

figure;plot(time(1:numFr,:),CumulDist(1:numFr,:),'-r') ; title ('Cumldispl vs. time');
hold on
axes('position', [.75 .15 .15 .15]);
plot(SpineData.SpineList(:,1), SpineData.SpineList(:,2), 'r');
saveas (gcf, [DataDir RecentFldr 'CumlPathTraveled' num2str(numFr)], 'pdf')




%plot neck movement
figure; plot(CurveMtx(3,1:numFr));
saveas (gcf, [DataDir RecentFldr 'NeckMovement' num2str(numFr)], 'pdf')
figure; plot(CurveMtx(3,:));
saveas (gcf, [DataDir RecentFldr 'NeckMovement'], 'pdf')


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