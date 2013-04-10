%% GET FINAL.mat file names if they exist
%G. kleemann
%eg, specify the suffixes for the folder and file
%[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'filt', 'final')

%in "folder" specify a current folder or 'recent' to pick the most recent
%RecentFldr will corespond the the folder name in specifed by "folder"
function [DateFileNms RecentFldr namedate] = GetTopFoldersMat (DataDir, folderID, fileID, folder)

fldrPath=[DataDir, filesep, 'RESULTS']
dirOutput = dir(fullfile(fldrPath, ['*RUN' folderID])); %specifiy source folder

%In case there is no "results" folder
if size(dirOutput, 1) < 1
fldrPath=DataDir
dirOutput = dir(fullfile([fldrPath,'*RUN' folderID])); %specifiy source folder
end 

if size(dirOutput, 1) < 1
fldrPath=DataDir
dirOutput = dir(fullfile([fldrPath,filesep, '*RUN' folderID])); %specifiy source folder
end 

%%
DateFldrNms = {dirOutput.name}';
if isempty(DateFldrNms); error (['I dont find your' folderID 'folder, make sure you built it']); end
namedate=sortrows([{dirOutput.datenum}', DateFldrNms]); %order by folder age

%% folder to work on
if strcmpi (folder, 'recent')
RecentFldr=namedate{end,2};%select the last folder built
else RecentFldr=folder{:,:};
end


%%
dirOutput = dir(fullfile([fldrPath, filesep, RecentFldr], ['*' fileID '.mat']))
DateFileNms = {dirOutput.name}';
end