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
[DateFileNms RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', 'recent')

for y=1:length(namedate(:,2));%cycle folders
    
    %% GET list of "final.mat" files in folder for looping
    [DateFileNms RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', namedate(y,2));
    
    %% stack the data from "final.mat" files
    %preallocate Matricies
    
    leng=length(DateFileNms)
    
    distanceMv=zeros(leng, 1);   %velocity=[velocity;vel];
    centroid=zeros(leng, 2);   %centroid=[centroid;CurrCent];
    imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
    CurveMtx=zeros(numpts-2, leng);%CurveMtx=[CurveMtx,SpineData.AngleLs];
    SpineList=cell(leng, 1);   %SpineList=[SpineList,SpineData.SpineList];
    Pointlist=cell(leng, 2);   %Pointlist=[Pointlist,SpineData.Pointlist];
    time=zeros(leng, 1);   %time=[time;insttime+time(end,:)];
    
%% check order, sort and add spacers for skipped or missing fra

[NameList]= GetImgNumOrdr(DateFileNms, '-', 'final')
 
for w=1:length(NameList);
        w
        NameList{w,1}
        load([DataDir filesep RecentFldr filesep DateFileNms{w,1}]);
        img1=varStruct.images.img1;
        [CurrCent]=FindCentr(varStruct.analysis.Img_Propfilt, 'CtrMass'); %CtrMass
        
        insttime=1/framerate; %seconds
        
        %use last centroid and current centroid to calcualte velocity
        if w < 2
            distanceMv(w,1)=0;
        else
            maxrow=max(find(centroid(:,1) >1));
            distMv=(pdist2(centroid(maxrow,:),CurrCent));%*MicroM_Pixel)/insttime;
            distanceMv(w,1)=distMv;
        end
        
        %cocatenate the order corrected matricies
        centroid(w,1:2)=CurrCent;
        imgs{w}=varStruct.images.imgBWL;%not sure if necc
        CurveMtx(:,w)=varStruct.SpineData.AngleLs;
        SpineList{w}=varStruct.SpineData.SpineList;
        Pointlist{w}=varStruct.SpineData.Pointlist;
        
        if w==1
            time(w)=0;
        else
            time(w)=insttime+time(w-1);
        end
        
        
        
        %% make figs
        
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
    
    %% SAVE EVERYTHING

    outputDir=[Alldata filesep 'ResultsFiles']
    mkdir(outputDir);
    NameOut=[outputDir filesep RecentFldr(1:end-8)];
    
    save ([NameOut 'SummData'], 'centroid', 'imgs', 'CurveMtx', 'SpineList', 'Pointlist', 'distanceMv', 'time', 'img1'); %does not save sub structure
    TrialData=dataset(centroid (:,1), centroid (:,2), distanceMv, time, 'VarNames', {'X','Y','Dist','time (sec)'});
    %export(TrialData, 'XLSfile', [DataDir RecentFldr 'movedata']);
    export(TrialData, 'XLSfile', [NameOut '_movedata']);
    
    CurveData=dataset(time, CurveMtx');
    %export(CurveData, 'XLSfile', [DataDir RecentFldr 'curvedata']);
    export(TrialData, 'XLSfile', [NameOut 'curvedata']);

%% 
    close all

    numFr =100; % frame count for the zoomed in view
    numFr2 =50;
    clims = [ -45 45 ];


    %full color full scale
    figure;imagesc(CurveMtx); colorbar
    saveas (gcf, [NameOut 'CurveMtxrfull'], 'pdf')
    
    MakeTricolMap
    %figure;imagesc(CurveMtx, clims)
    %colormap(RWBcMap2); colorbar
    %saveas (gcf, [NameOut 'CurveMtx'], 'pdf')
    
    figure;imagesc(CurveMtx(:,1:numFr), clims)
    colormap(RWBcMap2); colorbar
    saveas (gcf, [NameOut 'CurveMtx' num2str(numFr)], 'pdf')
    

    figure;imagesc(CurveMtx(:,1:numFr2), clims)
    colormap(RWBcMap2); colorbar
    saveas (gcf, [NameOut 'CurveMtx' num2str(numFr2)], 'pdf')
    
    figure;plot(centroid(:,1),centroid(:,2), '--rs',...
        'LineWidth', 2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g', 'MarkerSize', 1);...
        title ('centroid position'); xlim ([1,size(img1,2)]); ylim([1, size(img1, 1)]);
    saveas (gcf, [NameOut 'PathTraveled'], 'pdf')
    
%     figure;imagesc(CurveMtx)
%     figure;plot(centroid(1:numFr,1),centroid(1:numFr,2), '--rs',...
%         'LineWidth', 2, 'MarkerEdgeColor','k',...
%         'MarkerFaceColor','g', 'MarkerSize', 1);...
%         title ('centroid position'); xlim ([1,size(img1,2)]); ylim([1, size(img1, 1)]);
%     saveas (gcf, [NameOut 'PathTraveled' num2str(numFr)], 'pdf')
%     
%     
    
    %figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,450]); ylim([1,400]);
    
    %plot displacement
    figure;plot(time,distanceMv,'-r') ; title ('displacement vs. time');
    hold on
    axes('position', [.75 .15 .15 .15]);
    plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    saveas (gcf, [NameOut 'Distvstm'], 'pdf')
    
 %   figure;plot(time(1:numFr,:),distanceMv(1:numFr,:),'-r') ; title ('displacement vs. time');
 %   hold on
  %  axes('position', [.75 .15 .15 .15]);
 %   plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
  %  saveas (gcf, [NameOut 'DistVstm' num2str(numFr)], 'pdf')
    
    %Cummulative Dist
    CumulDist=cumsum(distanceMv)
    figure;plot(time,CumulDist,'-r') ; title ('Cumldispl vs. time');
    hold on
    axes('position', [.75 .15 .15 .15]);
    plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    saveas (gcf, [NameOut 'CumlPathTraveled'], 'pdf')
    
   % figure;plot(time(1:numFr,:),CumulDist(1:numFr,:),'-r') ; title ('Cumldispl vs. time');
   % hold on
   % axes('position', [.75 .15 .15 .15]);
   % plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
  %  saveas (gcf, [NameOut 'CumlPathTraveled' num2str(numFr)], 'pdf')

    %plot neck movement
    figure; plot(CurveMtx(3,1:numFr));
    saveas (gcf, [NameOut 'NeckMovement' num2str(numFr)], 'pdf')
    figure; plot(CurveMtx(3,:));
    saveas (gcf, [NameOut 'NeckMovement'], 'pdf')
    
    
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
    
    
 