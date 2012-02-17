%% GET FINAL.mat file names if they exist
%G. kleemann
%eg, specify the suffixes for the folder and file
%[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'filt', 'final')


function [DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, folderID, fileID)

fldrPath=[DataDir, filesep, 'RESULTS']
dirOutput = dir(fullfile(fldrPath, ['*RUN' folderID])); %specifiy source folder

%In case there is no "results" folder
if size(dirOutput, 1) < 1
fldrPath=DataDir
dirOutput = dir(fullfile([fldrPath,'*RUN' folderID])); %specifiy source folder
end 

DateFldrNms = {dirOutput.name}';
if isempty(DateFldrNms); error (['I dont find your' folderID 'folder, make sure you built it']); end
namedate=sortrows([{dirOutput.datenum}', DateFldrNms]); %order by folder age
RecentFldr=namedate{end,2};%select the last folder built
dirOutput = dir(fullfile([fldrPath, filesep, RecentFldr], ['*' fileID '.mat']))
DateFileNms = {dirOutput.name}';
end