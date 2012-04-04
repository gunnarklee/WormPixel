function ProcessSpine(varargin)

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
%OneWorm_CHR7Params


%% GET THE MOST RECENT "final.mat" folder for looping
[DateFileNms RecentFldr namedate] = GetFolderMat (DataDir, 'final', 'final');

%preallocate Matricies

for w=1:length(DateFileNms);
    %load([DataDir filesep 'RESULTS' filesep RecentFldr filesep DateFileNms{w,1}]);
    RUNfinalDir = [Alldata, 'RESULTS', filesep, TrialName, filesep, DateFldrNms{W} 'RUNfinal'];
    load([DataDir filesep RecentFldr filesep DateFileNms{w,1}]);

    % get spine angles
    varStruct.SpineData.AngleLs=GetAngles(varStruct.Pointlist.AngleLs)
    
saveThis([RUNfinalDir filesep SaveImNm, 'final.mat'], varStruct);%'ProcessDate'
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