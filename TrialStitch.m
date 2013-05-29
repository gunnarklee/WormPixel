function OneWormFigs(varargin)
%G. Kleemann - 5/28/13
%stitches together individual replicate (PIC_ folders) and sumarizes data by strain 
%   Detailed explanation goes here

%Version History

p = inputParser;
%p.addRequired('inputdir', @isdir);
p.addRequired('outputdir', @isdir);
p.addOptional('trialname', 'trialname', @ischar);

% Parse Inputs
try
    p.parse(varargin{:})
catch e1
    try
        %inputdir = uigetdir('choose the "Source" folder');
        outputdir = uigetdir('choose the "Destination" folder');
        p.parse(outputdir, varargin{:})
    catch e2
        exception = MException(...
            'twoimageDropCHR7_3Robot:arglist',...
            'Error in input argument list');
        exception = addCause(exception, e1);
        exception = addCause(exception, e2);
        throw(exception)
    end
end
%DataDir = p.Results.inputdir;
DataDir= [p.Results.outputdir filesep 'DoneRESULTS'];
Alldata = p.Results.outputdir;
TrialName = p.Results.trialname;
AlldataTop = Alldata(1:max(findstr(Alldata, filesep)))

%disp(['Input directory: ' DataDir]);
disp(['Output directory: ' Alldata]);
disp(['Trial name: ' TrialName]);

ProcessDate=(date);


%% Make new directories
% ErrorDir = [Alldata, 'RESULTS', filesep, TrialName, 'ErrorDir'];
% ProcessDate=(date);
% SCALETYPE='(Pixels)'; %unless calibrated to mm
% pixelPERmm=1; %default
% %SummaryFiles=[DataDir RecentFldr filesep 'SummData'],

%% SPECIFY PARAMETERS
OneWorm_CHR7Params;
% %% get PIXEL/MM CALIBRATION
% MeasureScale (resz, Alldata)

%% specifiy data matricies
scrsz = get(0,'ScreenSize');

Movement=[]; 
%Avg and variance %velocity, 
%aceleration, 
%Neck ampltude, 
%neck wavelength; 
%stroke/dist
SpineAngles=[]; 
%at static points on body 

%Wormimage=[]; % colect all worm images
%Worm Movie = image of worm translocating
%% update directories to find the results file
% [~, RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', 'recent');
% namedate

%load ([Alldata filesep 'Measure.mat']); %save parm -posctr s ellipse
%% get folder lists
%dirOutput = dir(fullfile([fldrPath, filesep, RecentFldr], ['*' fileID '.mat']))

%treatment names 
dirOutput = dir(fullfile([Alldata]))
TreatmentName = {dirOutput.name}'; 

%remove any directories with a dot ( these are not treament folders)
NoTrtm=strfind(TreatmentName, '.')
Index = find(cellfun('isempty', NoTrtm))
TreatmentName = TreatmentName(Index)

%% Cycle treatments
for y=1:length(TreatmentName);%cycle through treatments
    
%get replicate names and cycle through
dirOutput = dir(fullfile([Alldata, filesep, TreatmentName{y}, filesep, 'Done'], 'PIC_*'))
repName = {dirOutput.name}'; 

% cycle replicates concatenate and analyze data
for z=1:length(repName)
    dirOutput= dir(fullfile([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}]))
    repResult={dirOutput.name}';
    
    %hard code the output files > think of somthing more dynamic later
    SummData=load([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{11}])
     SwimStrokedata=importdata([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{12}])
    curvedata=importdata([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{13}])
    movedata=importdata([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{14}])

end
end

    %% GET list of "final.mat" files in folder for looping
    %[DateFileNms RecentFldr namedate] = GetTopFoldersMat (DataDir, 'final', 'final', namedate(y,2));
    
    %% check order, sort and add spacers for skipped or missing fra
    %Namelist has spacers, DateFileNmsDoes not
 %   try [NameList]= GetImgNumOrdr(DateFileNms, '-', 'final');
  %  catch
  %      [NameList]= GetImgNumOrdr(DateFileNms, '_', 'final');
   % end
    %% stack the data from "final.mat" files
    %preallocate Matricies
    leng=length(DateFileNms);  %UnPadded count
    PaddedLeng=length(NameList); %Padded count
    DateFileNms_Count=0%set up separate cont for UNPADEDDED LIST date file names
    insttime=1/framerate; %seconds
    
    distanceMv=zeros(PaddedLeng, 1);   %velocity=[velocity;vel];
    time=zeros(PaddedLeng, 1);   %time=[time;insttime+time(end,:)];
    centroid=zeros(PaddedLeng, 2);   %centroid=[centroid;CurrCent];
    %>>imgs=cell(leng, 1);  %imgs=[imgs;imgBWL];
    CurveMtx=zeros(numpts-2, 1);%CurveMtx=[CurveMtx,SpineData.AngleLs];
    SpineList=cell(PaddedLeng, 1);   %SpineList=[SpineList,SpineData.SpineList];
    Pointlist=cell(PaddedLeng, 2);   %Pointlist=[Pointlist,SpineData.Pointlist];
    
    for PaddedList_Count=1:length(NameList); %Put the data together
        PaddedList_Count
        NameList{PaddedList_Count,1}
        DateFileNms_Count=DateFileNms_Count+1;%set up separate cont for date file names
        fileDirect=[DataDir filesep RecentFldr]
        % check that here is a file (notspacer 'X') and it exists in a dir
        if and(isdir(fileDirect), ~strcmpi(NameList{PaddedList_Count,1}, 'X')); 
            load([fileDirect filesep NameList{PaddedList_Count,1}]);
            [CurrCent]=FindCentr(varStruct.analysis.Img_Propfilt, 'CtrMass'); %CtrMass
        end
        
        %% Add vales to the DATA matricies
        if NameList{PaddedList_Count,1} == 'X'; %then fill a blank in some matricies
            %keep zero tto retain spacerand goto next iteration
            %distanceMv(PaddedList_Count,1)=distanceMv(PaddedList_Count-1,1);
            time(PaddedList_Count) = insttime+max(time); %determined by FPS rate
            %centroid(PaddedList_Count,1:2)=[centroid(PaddedList_Count-1,1:2)]; %just take the last value
            continue;
        elseif PaddedList_Count == 1%Add actual value
            distanceMv(PaddedList_Count,1)=0;
            distMv=0;
            time(PaddedList_Count)=0; %determined by FPS rate
            centroid(PaddedList_Count,1:2)=CurrCent;
            %>imgs{DateFileNms_Count}=varStruct.images.imgBWL;%not sure if necc
            CurveMtx(:,PaddedList_Count)=varStruct.SpineData.AngleLs;
            SpineList{PaddedList_Count}=varStruct.SpineData.SpineList;
            Pointlist{PaddedList_Count}=varStruct.SpineData.Pointlist;
            %time(DateFileNms_Count)=insttime+time(DateFileNms_Count-1); %calculate times
            continue
        elseif (PaddedList_Count > 1);%after the first entry, put a VALUE or SPACER
            %CALCULATE VALUES TO BE PLOTTED
            %>f=varStruct.images.img1;
            %[CurrCent]=FindCentr(varStruct.analysis.Img_Propfilt, 'CtrMass'); %CtrMass
            %find the last centroid (highest ROW value)with a value - required since 0
            %values pad the unscored regions
            
            %% DISTANCE AND VELOCITY
            %use last centroid and current centroid to calcualte velocity
            % the odd loop is used because of discreopancies btwn padded and
            % unpadded list
            LastCentrRow=(max(find(centroid(:,1) >1)))-1
            %Find the last centroud psotion and distance moved
            if LastCentrRow > 1
                LastCentr=centroid(LastCentrRow-1,:);
                distMv=(pdist2(LastCentr,CurrCent))
            elseif isempty(LastCentrRow),
                distMv = 0;
            elseif LastCentrRow <= 1;
                distMv = 0;
            end
        end
 
    
    %end %add value or pad
    
    %% cocatenate the order corrected matricies
    %>imgs{DateFileNms_Count}=varStruct.images.imgBWL;%not sure if necc
    centroid(PaddedList_Count,1:2)=CurrCent;
    distanceMv(DateFileNms_Count,1)=distMv;
    CurveMtx(:,PaddedList_Count)=varStruct.SpineData.AngleLs;
    SpineList{PaddedList_Count}=varStruct.SpineData.SpineList;
    Pointlist{PaddedList_Count}=varStruct.SpineData.Pointlist;
    time(DateFileNms_Count)=insttime+time(DateFileNms_Count-1); %calculate times
    
    %%PLACE VUSIALIZE FUNCITON HERE
    %end % initialize (row 1) or pad
    %end % cycle through images - all images from one folder processed at this point
    
    end %image loop
    %% GET DERIVED PARAMETERS
    %these include path traveled (smoothed displacemnt)
    %and reset time for head subcycles (passed 0 degrees)
    
    % PATH PARAMETERS - USING ONLY NON ZERO VALUES 
    %Get the non zero values for centroid and corresonding time
    centroidPlot=centroid((find(~centroid(:,1)==0)), 1:2);
    
    % USE THIS FOR COLOR EVOLUTION
    timePlot=time(find(~centroid(:,1)==0));
    %[wkSpnX,wkspnY]=ind2sub(size(SpineData.SpineList), find(SpineData.SpineList));
    %cmap=colormap(jet (size(wkSpnX,1))) %copper
    
    %% smooth path with sliding window
    %STILL NEEDS WORK (run of ZEROS)
    % ** problem: path straight becomes static in the end (eg 37/265;% #223-258, jumps at 222 and 257)
  
    smoothWin=25;
    StraightWin=50;
    if length(centroidPlot)< smoothWin
        centroidPlot
        error(['usable data is shorter than the sliding window. Only ', num2str(length(centroidPlot)) ' datum in " ',namedate{y,2},' " folder. Refilter or add more pictures'])
    end    
    %Smooth forward
    %Last 25 (SmoothWin) are same value
    pathSmooth=[SmoothData(centroidPlot(:,1), smoothWin), SmoothData(centroidPlot(:,2), smoothWin)];
    %Last 50 (StraightWin) are same value
    pathStrgt=[SmoothData(centroidPlot(:,1), StraightWin), SmoothData(centroidPlot(:,2), StraightWin)];
    
    % % take last data point and draw a line to it with intermediate points
    % between
    % find the point where the movemb
    Breakpt=length(pathStrgt)-StraightWin*.1;
    pathStrgtStart=pathStrgt(1:Breakpt,1:2);
    %FILL IN  W/unsmoothed data rather than the single mean position (no
    %movement at all)
    pathStrgtEND=[];
    Startval=centroidPlot(Breakpt+1,1:2);
    Endval=centroidPlot(end,1:2);
    diffVal=Endval-Startval;
    intervals= length(centroidPlot)-(Breakpt+1);
    subdiv=diffVal/intervals;
    pathStrgtEND=[[Startval(1):subdiv(1):Endval(1)]', [Startval(2):subdiv(2):Endval(2)]'];
    
    %Stitch together
    pathStrgt=[pathStrgtStart;pathStrgtEND];
    
   %% VISUALIZE result
   % figure;plot(centroidPlot(2:end,1),centroidPlot(2:end,2), '-b',  'MarkerSize', 2);...
   %     hold on
    %plot(pathSmooth(:,1),pathSmooth(:,2), '--r',  'MarkerSize', 2);...
   % plot(pathStrgt(:,1),pathStrgt(:,2), '--k',  'MarkerSize', 2, 'lineWidth', 3);...
   %     pathSmooth(length(pathStrgt)-StraightWin+1:end,1)=centroidPlot(length(pathStrgt)-StraightWin+1:end, 1);
    %% TRAVEL Distance from STRAIGHTened path.
    %% problem: the smoothing window greates a dead space where movement is inappropeately clumped at the end of the trace
    % Solution :for now just chop off winsize portion ?check that speed calculations
    %are accurate still
    %slso could use actual displacement but this wil be inflated by
    %centroid waggle (not directed movement)
    
    TravelDist=zeros(size(centroid ,1),1);   %dont use pathstraight for size since it has NO GAPS;
   for WmPos=2:size(pathStrgt, 1); %Put the data together
        distmoved=pdist2(pathStrgt(WmPos-1,:), pathStrgt(WmPos,:));
        TravelDist(WmPos)=distmoved;
   end
    
    TravelDistTrim=TravelDist(1:(end-StraightWin-5)) % additional 5 is from top padding
    %% adjust to calibrated mm/pixel (1 if not calibrated)
    TravelDist = TravelDist/pixelPERmm;
    distanceMv = distanceMv/pixelPERmm;
    
    %% SWIM STROKE statistics
    ZeroLine=[zeros(length(time),1),time];
    OscilAngle=CurveMtx(3,1:end)';
    
    %STROKE COUNT
    %looking for crossing of the midline (angle changes from + to -
    %** has to be sustained on one side for more than 1 timepoint**
    [StrokeMk,  signlist] = SignChangeDetect(OscilAngle, time, insttime);
    
    StrokeCt=length(StrokeMk)-1;
    DistPerStroke=sum(TravelDist)/StrokeCt;
    
    %preallocate matricies
    InterStrokeTime=zeros(length(StrokeMk),1);
    AmplStroke=zeros(length(StrokeMk),1);
    
    for CurrStk=2:length(StrokeMk);
        StrokeStartTime=StrokeMk(CurrStk-1);
        StrokeEndTime=StrokeMk(CurrStk);
        StrokeStartROW=min(find(time >= StrokeMk(CurrStk-1)));
        StrokeEndROW=min(find(time >= StrokeMk(CurrStk)));
        %if thestroke end is not found find the time most similar to the stroke row end
        if isempty(StrokeEndROW);
            StrokeEndROW = MostSim(StrokeMk(CurrStk),time)
        end
        InterStrokeTime(CurrStk)=StrokeEndTime-StrokeStartTime;
        
        %Max angle calc
        if signlist(CurrStk) > 0
            %StrokePolarity(CurrStk)=1;
            AmplStroke(CurrStk)=max(OscilAngle(StrokeStartROW:StrokeEndROW));
        elseif signlist(CurrStk) < 0
            CurrStk
            StrokeStartROW
            StrokeEndROW
            size(OscilAngle)
            size(AmplStroke)
            AmplStroke(CurrStk)=min(OscilAngle(StrokeStartROW:StrokeEndROW));
            %StrokePolarity(CurrStk)=-1
            continue
        end
    end
    
    %% SAVE EVERYTHING
    outputDir=[Alldata filesep 'ResultsFiles'];
    mkdir(outputDir);;
    NameOut=[outputDir filesep RecentFldr(1:end-8) filesep];
    mkdir(NameOut);
    
    %MOVEMENT DATA
    save ([NameOut 'SummData'], 'centroid', 'imgs', 'CurveMtx',...
        'SpineList', 'Pointlist', 'distanceMv', 'time',... % 'img1',
        'InterStrokeTime', 'DistPerStroke',...
        'AmplStroke', 'StrokeMk', 'signlist') %does not save sub structure
    %TrialData=dataset(centroid (:,1), centroid (:,2), distanceMv, time,..
    %'VarNames', {'X','Y','Dist','time (sec)'});
    TrialData=dataset(centroid (:,1), centroid (:,2), distanceMv, time,...
        'VarNames', {'X','Y','Dist','time (sec)'});
    %export(TrialData, 'XLSfile', [DataDir RecentFldr 'movedata']);
  %  if ispc
   %     export(TrialData, 'XLSfile', [NameOut 'movedata']);
   % else
        export(TrialData,'file', [NameOut 'movedata'],'Delimiter',',');
    %end
    
    %STROKE DATA
    Strokedata=dataset(InterStrokeTime, AmplStroke, StrokeMk,...
        signlist, 'VarNames', {'StrokeDurn', 'StrokeAmpl',...
        'StrokeMk', 'signlist'}); %maybe add stroke specific power...'Strokeower'
    
    %export(TrialData, 'XLSfile', [DataDir RecentFldr 'movedata']);
   %if ispc
    %    export(Strokedata, 'XLSfile', [NameOut 'SwimStrokedata']);
   % else
        export(Strokedata,'file', [NameOut 'SwimStrokedata'],'Delimiter',',');
    %end
    
    %CURVE DATA
    CurveData=dataset(time, CurveMtx');
    %export(CurveData, 'XLSfile', [DataDir RecentFldr 'curvedata']);
    %if ispc
     %   export(CurveData, 'XLSfile', [NameOut 'curvedata']);
   % else
        export(CurveData,'file', [NameOut 'curvedata'],'Delimiter',',');
    %end
    
    %% PLOT FIGURES
    close all
    
    numFr =50; % frame count for the zoomed in view
    %>numFr2 =25;
    %>clims = [-35 35];
    
    %% curve matricies
    %full color full scale
    figure;imagesc(CurveMtx); colorbar
    saveas (gcf, [NameOut 'CurveMtxrfull'], 'pdf')
    
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
    
    %% Path traveled - COLOR EVOLVING PATH!
    cmapTmRG=[];
    cmap=colormap(jet (size(centroidPlot,1))); %copper
    cmapTime=colormap(jet (size(time,1))); %copper
    cmapTmRGB=cat(3,cmapTime(:,1), cmapTime(:,2), cmapTime(:,3)); % for time scale
    
   cmapTmRGB =repmat(cmapTmRGB, length(cmapTmRGB), 10); 
    
  % figure; imshow(cmapTmRGB);
%% 
    figure; plot(centroidPlot(1,1),centroidPlot(1,2), '-b',  'MarkerSize', 1);...
    title ('centroid position');    
    xlim ([1,4500]); ylim([1,2500]);
    hold on

    for Pt=2:length(centroidPlot)
    p=plot(centroidPlot(Pt-1:Pt,1),centroidPlot(Pt-1:Pt,2));... 
    set(p,'Color',cmap(Pt,:),'LineWidth',1)
    end
    
    %****REINSTATE***xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
    
    saveas (gcf, [NameOut 'Path of centroid across frames'], 'pdf')
    
    %% PLOT Smooth path -
    figure; plot(centroidPlot(1,1),centroidPlot(1,2), '-b',  'MarkerSize', 2);...
    title ('centroid position');    
    hold on

    for Pt=2:length(centroidPlot)
    p=plot(centroidPlot(Pt-1:Pt,1),centroidPlot(Pt-1:Pt,2));... 
    set(p,'Color',cmap(Pt,:),'LineWidth',1)
    end    
    
    %figure;plot(centroidPlot(2:end,1),centroidPlot(2:end,2), '-b',  'MarkerSize', 2);...
    %    hold on
    %plot(pathSmooth(:,1),pathSmooth(:,2), '--r',  'MarkerSize', 2);...
    plot(pathStrgt(:,1),pathStrgt(:,2), '--k',  'MarkerSize', 2, 'lineWidth', 2);...
        legend ('centroid path', 'path traveled (smoothed)')
    title ('travel path');
    %xlim ([1,size(img1,1)*(resz*1.25)]); ylim([1, size(img1,2)*(resz*1.25)]);
    saveas (gcf, [NameOut 'PathTraveled'], 'pdf')
    
    
    %% plot displacement
    
    % outrageously high values (user moved plate) are smothehd over
    % usin> 10 standard deviations
    distanceMv(find(distanceMv > (std(distanceMv)*10)))=mean(distanceMv)
    figure;plot(time(1:end-1),distanceMv(1:end-1),'-r') ; title ('Centroid displacement between frames');
    %last frame is zero,
    hold on
    xlabel ('time (Sec)');
    ylabel (['centorid movement', SCALETYPE]);
    %>>axes('position', [.75 .15 .15 .15]);
    %>>plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    %>>saveas (gcf, [NameOut 'centroid movement'], 'pdf')
    
    
    %% compensate for smoothing losses in path for plots
    %% problem: the smoothing window greates a dead space where movement is inappropeately clumped at the end of the trace
    % Solution :for now just chop off winsize portion ?check that speed calculations
    %are accurate still
    %could use actual displacement but this wil be inflated by centroid waggle (not directed movement)
   
    TravelDistTrim=TravelDist(1:(end-StraightWin-5)) % additional 5 is from top padding
    timeTrim=time(1:(end-StraightWin-5))
    
    %% plot PATH TRAVELED  displacement
    figure;plot(timeTrim(1:end), TravelDistTrim(1:end),'-r') ; title ('displacement between frames');
    %last frame is zero,
    %hold on
    xlabel ('time (Sec)');
    ylabel (['centorid movement', SCALETYPE]);
    %>>axes('position', [.75 .15 .15 .15]);
    %>>plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    %saveas (gcf, [NameOut 'DIstance Traveled'], 'pdf')
    
    %% Cummulative Dist
    CumulDist=cumsum(TravelDistTrim);
    figure;plot(timeTrim(1:length(CumulDist)),CumulDist,'-r');
    title (['Cumlative travel (straight path)' SCALETYPE ' vs. time']);
    xlabel ('time (Sec)');
    ylabel (['Cumulative Travel distance ', SCALETYPE]);
    
    %>> hold on
    %>> axes('position', [.75 .15 .15 .15]);
    %>> plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    saveas (gcf, [NameOut 'Cuml_Travel'], 'pdf')
    
    %% plot SUBSET neck movement
    
    SubsetStroke=StrokeMk(find(StrokeMk< max(time(1:numFr))))

    figure;
    plot(time(1:numFr),OscilAngle(1:numFr),'-b');
    hold on
    plot(ZeroLine((1:numFr),2),ZeroLine((1:numFr),1),'-k')
%%>>>lengths diff here somtimes!
    plot(SubsetStroke, 0, '*r')
        h=legend ('head angle','stroke reset line', 'stroke marker');
    set (h, 'Location', 'NorthEast');
    
    title (['Head oscillations' num2str(numFr) 'frames'])
    xlabel ('time (Sec)');
    ylabel ('degrees from the midline');
    
    saveas (gcf, [NameOut 'Neck Movement' num2str(numFr)], 'pdf');
    %% plot FULL RUN neck movement
    fullRun=length(time)-1; %puts zeros at tend
    
    figure;
    plot(time(1:fullRun),OscilAngle(1:fullRun),'-b');
    hold on
    plot(ZeroLine((1:fullRun),2),ZeroLine((1:fullRun),1),'-k')
    plot(StrokeMk, 0, '*r');
    h=legend ('head angle','stroke reset line', 'stroke marker');
    set (h, 'Location', 'SouthEast');
    
    
    title ('Head oscillations Full Run')
    xlabel ('time (Sec)');
    ylabel (['degrees from the midline']);
    xlim([0,max(time)]);
    
    saveas (gcf, [NameOut 'Neck MovementFullRun'], 'pdf')
    %% PLOT CYCLE (ENDS)
end %function 
  
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

function SimRow = MostSim(goal, targetMtx)
goalMtx=repmat(goal,size(targetMtx));
diffMtx=goalMtx-targetMtx;
SimRow=find(diffMtx==min(diffMtx));
end