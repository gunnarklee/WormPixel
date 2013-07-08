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

%observation matricies
    Trtm={};
    TrtmLong={};
    RepLong={};
    
%Single rep matricies 
%- mean standard deviation and absolue value mean(avg dist from 0) across single reps

    MoveDtMEAN_rep= [];
    MoveDtSTD_rep= []; 
    
    StrokeDtMEAN_rep= [];
    StrokeDtSTD_rep= [];
    
    CurveDtMEAN_rep= [];    
    CurveDtSTD_rep= [];
    CurveDtABSL_MEAN_rep=[];

%average for strain "overall" matricies.    
    MoveGrandMeans= [];
    MoveStdofmeans= []; %mean variance across replicates within each strain 
    MoveMeanofStd= []; %mean variance in replicates across strains
    
    StrokeGrandMeans= [];
    StrokeSTDofMeans= []; %mean variance across replicates within each strain 
    StrokeMeanofSTD= []; %mean variance in replicates across strains
    
    CurveGrandMeans= [];    
    CurveSTDofMeans= []; %mean variance across replicates within each strain  
    CurveMeanofSTD= [];  %mean variance in replicates across strains
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
 TreatmentName{y}   
%get replicate names and cycle through
dirOutput = dir(fullfile([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles'], 'PIC_*'));
repName = {dirOutput.name}'; 

% cycle replicates concatenate and analyze data
%these variables are reset for each strain
    MoveDtMEANtmp= [];
    MoveDtSTDtmp= [];
    
    StrokeDtMEANtmp= [];
    StrokeDtSTDtmp= [];
    
    CurveDtMEANtmp= [];    
    CurveDtSTDtmp= [];
    CurveDtABSL_MEANtmp=[];

    if size (repName, 1) == 0 ; continue; end 
for z=1:length(repName)
   %if size (repName, 1) == 0 ; continue; end 
    display (['working on rep' repName{z}, 'the', num2str(z), 'of', num2str(length(repName))]);
    dirOutput= dir(fullfile([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}]));
    repResult={dirOutput.name}';
    
% check for the required files
%% hard coded the output files > think of something more dynamic later
   CurrRep = [Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}]
   try
    SummData=load([CurrRep, filesep, repResult{11}]);
    catch
        error (['I can not find the', repResult{11} 'file, is it missing?'])
    end
    
    try
    SwimStrokedata = dataset('file',[CurrRep, filesep, repResult{12}],'delimiter',',','ReadVarNames',true);    
    %SwimStrokedata=importdata();
    %SwimStrokedata=dataset(importdata([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{12}]));
    catch
        error (['I can not find the', repResult{12} 'file, is it missing?'])
    end
    
    try
    curvedata = dataset('file',[CurrRep, filesep, repResult{13}],'delimiter',',','ReadVarNames',true);    
    %curvedata=importdata([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{13}]);
    catch
        error (['I can not find the', repResult{13} 'file, is it missing?'])
    end
    
    try
    movedata = dataset('file',[CurrRep, filesep, repResult{14}],'delimiter',',','ReadVarNames',true);      
    %movedata=importdata([Alldata, filesep, TreatmentName{y}, filesep, 'ResultsFiles', filesep, repName{z}, filesep, repResult{14}]);
    catch
        error (['I can not find the', repResult{14} 'file, is it missing?'])
    end
    
% Compile data into structures
%     curveDataSum.Headers.(TreatmentName{y}).(repName{z})=curvedata.colheaders;
%     curveDataSum.Mean.(TreatmentName{y}).(repName{z})=mean(curvedata.data(6:end, :));
%     curveDataSum.Min.(TreatmentName{y}).(repName{z})=min(curvedata.data(6:end, :));
%     curveDataSum.Max.(TreatmentName{y}).(repName{z})=max(curvedata.data(6:end, :));
%     curveDataSum.Range.(TreatmentName{y}).(repName{z})=range(curvedata.data(6:end, :));
%     curveDataSum.StdDev.(TreatmentName{y}).(repName{z})=std(curvedata.data(6:end, :));
%     curveDataSum.total.(TreatmentName{y}).(repName{z})=sum(curvedata.data(6:end, :));
%     
    SummDataSum.(TreatmentName{y}).(repName{z})=SummData;

%     SwimStrokeSum.Headers.(TreatmentName{y}).(repName{z})=SwimStrokedata.colheaders;
%     SwimStrokeSum.Mean.(TreatmentName{y}).(repName{z})=mean(SwimStrokedata.data(6:end, :));
%     SwimStrokeSum.Min.(TreatmentName{y}).(repName{z})=min(SwimStrokedata.data(6:end, :));
%     SwimStrokeSum.Max.(TreatmentName{y}).(repName{z})=max(SwimStrokedata.data(6:end, :));
%     SwimStrokeSum.Range.(TreatmentName{y}).(repName{z})=range(SwimStrokedata.data(6:end, :));
%     SwimStrokeSum.StdDev.(TreatmentName{y}).(repName{z})=std(SwimStrokedata.data(6:end, :));
%     SwimStrokeSum.total.(TreatmentName{y}).(repName{z})=sum(SwimStrokedata.data(6:end, :));
%    
%     movedataSum.Headers.(TreatmentName{y}).(repName{z})=movedata.colheaders;
%     movedataSum.Mean.(TreatmentName{y}).(repName{z})=mean(movedata.data(6:end, :));
%     movedataSum.Min.(TreatmentName{y}).(repName{z})=min(movedata.data(6:end, :));
%     movedataSum.Max.(TreatmentName{y}).(repName{z})=max(movedata.data(6:end, :));
%     movedataSum.Range.(TreatmentName{y}).(repName{z})=range(movedata.data(6:end, :));
%     movedataSum.StdDev.(TreatmentName{y}).(repName{z})=std(movedata.data(6:end, :));
%     movedataSum.total.(TreatmentName{y}).(repName{z})=sum(movedata.data(6:end, :));
%     
%% temporary strain specific replicates files
    
    MoveDtMEANtmp= [MoveDtMEANtmp; mean(double(movedata(6:end,:)))];
    MoveDtSTDtmp= [MoveDtSTDtmp; std(double(movedata(6:end, :)))];
    
    StrokeDtMEANtmp= [StrokeDtMEANtmp; mean(double(SwimStrokedata(6:end, :)))];
    StrokeDtSTDtmp= [StrokeDtSTDtmp; std(double(SwimStrokedata(6:end, :)))];
    
    CurveDtMEANtmp= [CurveDtMEANtmp; mean(double(curvedata(6:end, :)))]    ;
    CurveDtABSL_MEANtmp= [CurveDtABSL_MEANtmp; mean(abs(double(curvedata(6:end, :))))] 
    CurveDtSTDtmp= [CurveDtSTDtmp; std(double(curvedata(6:end, :)))];
    
    TrtmLong=[TrtmLong; TreatmentName{y}];
    RepLong=[RepLong; repName{z}];
   % export(CurveData,'file', [NameOut 'curvedata'],'Delimiter',',')
end

%% COMPILE data into matricies

% compile rep data
    MoveDtMEAN_rep= [MoveDtMEAN_rep; MoveDtMEANtmp];
    MoveDtSTD_rep= [MoveDtSTD_rep;MoveDtSTDtmp];
    
    StrokeDtMEAN_rep= [StrokeDtMEAN_rep;StrokeDtMEANtmp];
    StrokeDtSTD_rep= [StrokeDtSTD_rep;StrokeDtSTDtmp];
    
    CurveDtMEAN_rep= [CurveDtMEAN_rep;CurveDtMEANtmp]    ;
    CurveDtSTD_rep= [CurveDtSTD_rep;CurveDtSTDtmp];
    CurveDtABSL_MEAN_rep=[CurveDtABSL_MEAN_rep;CurveDtABSL_MEANtmp]
        
% Average rep data to get strain data 
    Trtm=[Trtm; TreatmentName{y}]
    
    MoveGrandMeans= [MoveGrandMeans; mean(MoveDtMEANtmp)];
    MoveStdofmeans= [MoveStdofmeans; std(MoveDtMEANtmp)];
    MoveMeanofStd= [MoveMeanofStd; mean(MoveDtSTDtmp)];
    
    StrokeGrandMeans= [StrokeGrandMeans; mean(StrokeDtMEANtmp)];
    StrokeSTDofMeans= [StrokeSTDofMeans; std(StrokeDtMEANtmp)];
    StrokeMeanofSTD= [StrokeMeanofSTD; mean(StrokeDtSTDtmp)]   ;
    
    CurveGrandMeans= [CurveGrandMeans; mean(CurveDtMEANtmp)]    ;
    CurveSTDofMeans= [CurveSTDofMeans; std(CurveDtMEANtmp)];
    CurveMeanofSTD= [CurveMeanofSTD; mean(CurveDtSTDtmp)];
       
end
%% get headers
    Moveheaders=movedata.Properties.VarNames
    Strokeheaders=SwimStrokedata.Properties.VarNames
    Curveheaders=curvedata.Properties.VarNames
    
%% Make into DATASETS  with headers and obs names
    
    
    MoveDtMEAN_rep = [dataset(TrtmLong), dataset(RepLong), mat2dataset(MoveDtMEAN_rep, 'VarNames', Moveheaders)];
    MoveDtSTD_rep = [dataset(TrtmLong), dataset(RepLong), mat2dataset(MoveDtSTD_rep, 'VarNames', Moveheaders)];
    
    StrokeDtMEAN_rep = [dataset(TrtmLong), dataset(RepLong), mat2dataset(StrokeDtMEAN_rep , 'VarNames', Strokeheaders)];
    StrokeDtSTD_rep = [dataset(TrtmLong), dataset(RepLong), mat2dataset(StrokeDtSTD_rep, 'VarNames', Strokeheaders)];
    
    CurveDtMEAN_rep = [dataset(TrtmLong), dataset(RepLong), mat2dataset(CurveDtMEAN_rep(:,2:end), 'VarNames', Curveheaders(2:end))];
    CurveDtSTD_rep = [dataset(TrtmLong), dataset(RepLong), mat2dataset(CurveDtSTD_rep(:,2:end) , 'VarNames', Curveheaders(2:end))];
    CurveDtABSL_MEAN_rep= [dataset(TrtmLong), dataset(RepLong), mat2dataset(CurveDtABSL_MEAN_rep(:,2:end), 'VarNames', Curveheaders(2:end))];
    
    %     MoveMeans_Overall=mat2dataset(MoveGrandMeans, 'VarNames', Moveheaders, 'ObsNames', Trtm);
%     MoveStdofmeans_Overall= mat2dataset(MoveStdofmeans, 'VarNames', Moveheaders, 'ObsNames', Trtm);
%     MoveMeanofStd_Overall= mat2dataset(MoveMeanofStd, 'VarNames', Moveheaders, 'ObsNames', Trtm);
% 
%     StrokeMeans_Overall = mat2dataset(StrokeGrandMeans, 'VarNames', Strokeheaders, 'ObsNames', Trtm);
%     StrokeSTDofMeans_Overall = mat2dataset(StrokeSTDofMeans, 'VarNames', Strokeheaders, 'ObsNames', Trtm);
%     StrokeMeanofSTD_Overall = mat2dataset(StrokeMeanofSTD, 'VarNames', Strokeheaders, 'ObsNames', Trtm)  ;
%     
%     CurveMeans_Overall = mat2dataset(CurveGrandMeans, 'VarNames', Curveheaders, 'ObsNames', Trtm);
%     CurveSTDofMeans_Overall =  mat2dataset(CurveSTDofMeans, 'VarNames', Curveheaders, 'ObsNames', Trtm);
%     CurveMeanofSTD_Overall =  mat2dataset(CurveMeanofSTD, 'VarNames', Curveheaders, 'ObsNames', Trtm)


%% Make overall summaries into datasets with headers and obs names
    MoveMeans_Overall=mat2dataset(MoveGrandMeans, 'VarNames', Moveheaders, 'ObsNames', Trtm)
    MoveStdofmeans_Overall= mat2dataset(MoveStdofmeans, 'VarNames', Moveheaders, 'ObsNames', Trtm);
    MoveMeanofStd_Overall= mat2dataset(MoveMeanofStd, 'VarNames', Moveheaders, 'ObsNames', Trtm);

    StrokeMeans_Overall = mat2dataset(StrokeGrandMeans, 'VarNames', Strokeheaders, 'ObsNames', Trtm);
    StrokeSTDofMeans_Overall = mat2dataset(StrokeSTDofMeans, 'VarNames', Strokeheaders, 'ObsNames', Trtm)
    StrokeMeanofSTD_Overall = mat2dataset(StrokeMeanofSTD, 'VarNames', Strokeheaders, 'ObsNames', Trtm)  
    
    CurveMeans_Overall = mat2dataset(CurveGrandMeans(:,2:end), 'VarNames', Curveheaders(2:end), 'ObsNames', Trtm);
    CurveSTDofMeans_Overall =  mat2dataset(CurveSTDofMeans(:,2:end), 'VarNames', Curveheaders(2:end), 'ObsNames', Trtm)
    CurveMeanofSTD_Overall =  mat2dataset(CurveMeanofSTD(:,2:end), 'VarNames', Curveheaders(2:end), 'ObsNames', Trtm)

%% MAKE OUTPUT TEXT FILES
    % compiled data, CSV/dataset file of summary data
    outputDir=[Alldata filesep 'ResultsFiles'];
    mkdir(outputDir);
    NameOut=[outputDir filesep];
    mkdir(NameOut);

    %export data sheets
    %export(TrtmtSum,'file', [NameOut 'movedata'],'Delimiter',',')
    export(MoveDtMEAN_rep,'file', [NameOut 'MoveDtMEAN_rep'],'Delimiter',',');
    %export(MoveDtSTD_rep,'file', [NameOut 'MoveDtSTD_rep'],'Delimiter',',');
    export(StrokeDtMEAN_rep,'file', [NameOut 'StrokeDtMEAN_rep'],'Delimiter',',');
    %export(StrokeDtMEAN_rep,'file', [NameOut 'StrokeDtMEAN_rep'],'Delimiter',',');
    export(CurveDtMEAN_rep,'file', [NameOut 'CurveDtMEAN_rep'],'Delimiter',',');
    export(CurveDtSTD_rep,'file', [NameOut 'CurveDtSTD_rep'],'Delimiter',',');
    export(MoveMeans_Overall,'file', [NameOut 'MoveMeans_Overall'],'Delimiter',',');
    %export(MoveStdofmeans_Overall,'file', [NameOut 'curvedata'],'Delimiter',',');
    export(StrokeMeans_Overall,'file', [NameOut 'StrokeMeans_Overall'],'Delimiter',',');
    %export(StrokeSTDofMeans_Overall,'file', [NameOut 'StrokeSTDofMeans_Overall'],'Delimiter',',');
    export(CurveMeans_Overall,'file', [NameOut 'CurveMeans_Overall'],'Delimiter',',');
    export(CurveSTDofMeans_Overall,'file', [NameOut 'CurveSTDofMeans_Overall'],'Delimiter',',');
    %export(CurveData,'file', [NameOut 'curvedata'],'Delimiter',',');


%% GRAPH IT - summaries into datasets with headers and obs names

% Plot average angles along body
%get movement data for indivuduals split by strian - not sure about these
%units

% PLOT Average angles Average body position
varNames=CurveMeans_Overall.Properties.VarNames(1:end)'
Dt=[CurveMeans_Overall.Var2_1, CurveMeans_Overall.Var2_2, CurveMeans_Overall.Var2_3,...
CurveMeans_Overall.Var2_4, CurveMeans_Overall.Var2_5, CurveMeans_Overall.Var2_6,...
CurveMeans_Overall.Var2_7, CurveMeans_Overall.Var2_8, CurveMeans_Overall.Var2_9,...
CurveMeans_Overall.Var2_10, CurveMeans_Overall.Var2_11];
Y=CurveMeans_Overall.Properties.ObsNames

%specify colors
set(gca,'NextPlot','replacechildren');
set(gca,'ColorOrder', lines(size(Y,1)));
figure; parallelcoords(Dt, 'group',Y, 'labels', varNames);
% Legend replot
[legend_h,object_h,plot_h,text_strings] = legend()
legend(plot_h, text_strings, 'Location', 'best')

saveas (gcf, [NameOut 'CurveAvgs_overall'], 'pdf')
    

% gplotmatrix(ratings(:,1:2),ratings(:,[4 7]),group,... 
%             'br','.o',[],'on','',categories(1:2,:),... 
%              categories([4 7],:))


%% Plot body angle for REPs
varNames=CurveDtMEAN_rep.Properties.VarNames(3:end)'
Dt=[CurveDtMEAN_rep.Var2_1, CurveDtMEAN_rep.Var2_2, CurveDtMEAN_rep.Var2_3,...
CurveDtMEAN_rep.Var2_4, CurveDtMEAN_rep.Var2_5, CurveDtMEAN_rep.Var2_6,...
CurveDtMEAN_rep.Var2_7, CurveDtMEAN_rep.Var2_8, CurveDtMEAN_rep.Var2_9,...
CurveDtMEAN_rep.Var2_10, CurveDtMEAN_rep.Var2_11];
Y=CurveDtMEAN_rep.TrtmLong

%specify colors
set(gca,'NextPlot','replacechildren');
set(gca,'ColorOrder', lines(size(Y,1)));
figure; parallelcoords(Dt, 'group',Y, 'labels', varNames);
% Legend replot
[legend_h,object_h,plot_h,text_strings] = legend()
legend(plot_h, text_strings, 'Location', 'best')

saveas (gcf, [NameOut 'CurveAvgs_reps'], 'pdf')
    
%figure;andrewsplot(Dt, 'group',Y);
%figure;hbar(Dt, 'group',Y, 'labels', varNames);

MakeTricolMap
colormap(RWBcMap2)


%% plot as heatmap         
figure;imagesc(double(CurveMeans_Overall(:,1:end))); colorbar %specifiy title and axis
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
set(gca,'YTick',[1:length(trmt)])
set(gca,'YTickLabel', trmt)
saveas (gcf, [NameOut 'CurveAvgsHeat_Overall'], 'pdf')


figure;imagesc(double(CurveDtMEAN_rep(:,3:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
YtickInt=[1:RepCt/length(trmt):RepCt]
set(gca,'YTick',[1:YtickInt])
set(gca,'YTickLabel', trmt)
saveas (gcf, [NameOut 'CurveAvgsHeat_reps'], 'pdf')


figure;imagesc(double(CurveDtMEAN_rep(:,3:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
YtickInt=[1:RepCt/length(trmt):RepCt]
set(gca,'YTick',[1:YtickInt])
set(gca,'YTickLabel', trmt)
saveas (gcf, [NameOut 'CurveAvgsHeat_reps'], 'pdf')


figure;imagesc(double(CurveDtABSL_MEAN_rep(:,3:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
YtickInt=[1:RepCt/length(trmt):RepCt]
set(gca,'YTick',[1:YtickInt])
set(gca,'YTickLabel', trmt)
saveas (gcf, [NameOut 'CurveDtABSL_MEAN_rep'], 'pdf')


figure;imagesc(double(CurveSTDofMeans_Overall(:,1:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
YtickInt=[1:RepCt/length(trmt):RepCt]
set(gca,'YTick',[1:YtickInt])
set(gca,'YTickLabel', trmt)
saveas (gcf, [NameOut 'CurveSTDofMeans_Overall'], 'pdf')


figure;imagesc(double(CurveMeanofSTD_Overall(:,1:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
YtickInt=[1:RepCt/length(trmt):RepCt]
set(gca,'YTick',[1:YtickInt])
set(gca,'YTickLabel', trmt)

saveas (gcf, [NameOut 'CurveMeanofSTD_Overall'], 'pdf')


figure;imagesc(double(CurveMeanofSTD_Overall(:,1:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
set(gca,'YTick',[1:length(trmt)])
set(gca,'YTickLabel', trmt)
title('avg variance(std) per strain')
saveas (gcf, [NameOut 'CurveMeanofSTD_Overall'], 'pdf')

figure;imagesc(double(CurveDtSTD_rep(:,3:end))); colorbar %specifiy title and axis
colormap(RWBcMap2)
trmtTbl=tabulate(CurveDtMEAN_rep.TrtmLong)
trmt=trmtTbl(:,1)
RepCt=size(CurveDtMEAN_rep, 1)
YtickInt=[1:RepCt/length(trmt):RepCt]
set(gca,'YTick',[1:YtickInt])
set(gca,'YTickLabel', trmt)
saveas (gcf, [NameOut 'CurveDtSTD_rep'], 'pdf')
%% plot and save summary output

%montage individual curve figures

        %% Add vales to the DATA matricies
%         if NameList{PaddedList_Count,1} == 'X'; %then fill a blank in some matricies
%             %keep zero tto retain spacerand goto next iteration
%             %distanceMv(PaddedList_Count,1)=distanceMv(PaddedList_Count-1,1);
%             time(PaddedList_Count) = insttime+max(time); %determined by FPS rate
%             %centroid(PaddedList_Count,1:2)=[centroid(PaddedList_Count-1,1:2)]; %just take the last value
%          
%         elseif PaddedList_Count == 1%Add actual value
%             distanceMv(PaddedList_Count,1)=0;
%             distMv=0;
%             time(PaddedList_Count)=0; %determined by FPS rate
%             centroid(PaddedList_Count,1:2)=CurrCent;
%             %>imgs{DateFileNms_Count}=varStruct.images.imgBWL;%not sure if necc
%             CurveMtx(:,PaddedList_Count)=varStruct.SpineData.AngleLs;
%             SpineList{PaddedList_Count}=varStruct.SpineData.SpineList;
%             Pointlist{PaddedList_Count}=varStruct.SpineData.Pointlist;
%             %time(DateFileNms_Count)=insttime+time(DateFileNms_Count-1); %calculate times
%             continue
%         elseif (PaddedList_Count > 1);%after the first entry, put a VALUE or SPACER
            %CALCULATE VALUES TO BE PLOTTED
            %>f=varStruct.images.img1;
            %[CurrCent]=FindCentr(varStruct.analysis.Img_Propfilt, 'CtrMass'); %CtrMass
            %find the last centroid (highest ROW value)with a value - required since 0
            %values pad the unscored regions
            
%             %% DISTANCE AND VELOCITY
%             %use last centroid and current centroid to calcualte velocity
%             % the odd loop is used because of discreopancies btwn padded and
%             % unpadded list
%             LastCentrRow=(max(find(centroid(:,1) >1)))-1
%             %Find the last centroud psotion and distance moved
%             if LastCentrRow > 1
%                 LastCentr=centroid(LastCentrRow-1,:);
%                 distMv=(pdist2(LastCentr,CurrCent))
%             elseif isempty(LastCentrRow),
%                 distMv = 0;
%             elseif LastCentrRow <= 1;
%                 distMv = 0;
%             end
%         end
 
    
    %end %add value or pad
    
    %% cocatenate the order corrected matricies
    %>imgs{DateFileNms_Count}=varStruct.images.imgBWL;%not sure if necc
%     centroid(PaddedList_Count,1:2)=CurrCent;
%     distanceMv(DateFileNms_Count,1)=distMv;
%     CurveMtx(:,PaddedList_Count)=varStruct.SpineData.AngleLs;
%     SpineList{PaddedList_Count}=varStruct.SpineData.SpineList;
%     Pointlist{PaddedList_Count}=varStruct.SpineData.Pointlist;
%     time(DateFileNms_Count)=insttime+time(DateFileNms_Count-1); %calculate times
%     
    %%PLACE VUSIALIZE FUNCITON HERE
    %end % initialize (row 1) or pad
    %end % cycle through images - all images from one folder processed at this point
  
    %% GET DERIVED PARAMETERS
    %these include path traveled (smoothed displacemnt)
    %and reset time for head subcycles (passed 0 degrees)
    
    % PATH PARAMETERS - USING ONLY NON ZERO VALUES 
    %Get the non zero values for centroid and corresonding time
%     centroidPlot=centroid((find(~centroid(:,1)==0)), 1:2);
%     
%     % USE THIS FOR COLOR EVOLUTION
%     timePlot=time(find(~centroid(:,1)==0));
%     %[wkSpnX,wkspnY]=ind2sub(size(SpineData.SpineList), find(SpineData.SpineList));
%     %cmap=colormap(jet (size(wkSpnX,1))) %copper
    
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
    mkdir(outputDir);
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