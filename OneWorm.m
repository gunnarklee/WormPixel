function OneWorm(varargin)
%G. Kleemannn 5/2013
%OneWorm.m - start point for One Worm analysis pipeline

% Parser for case when inputdir and outputdir are specified
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
Alldata = p.Results.outputdir;
TrialName = p.Results.trialname;
AlldataTop = Alldata(1:max(findstr(Alldata, filesep)))

%disp(['Input directory: ' DataDir]);
disp(['Output directory: ' Alldata]);
disp(['Trial name: ' TrialName]);

%% WORKFLOW for OneWorm.
% specify your path to the input directory. 
%InputDir='C:\Documents and Settings\mlab\Desktop\reproductive_span_robot\RSbatch2screen1-verify12_2011'
%..InputDir='/Users/budoday/Dropbox/Movement/Movement_Analysis/REDO'

% On a PC always specify a trialName

%InputDir='C:\Documents and Settings\mlab\My Documents\Dropbox\WormMovement\AdulanteSwimming\SwimData\fld2'

%Documents\Dropbox\WormMovement\AdulanteSwimming\SwimData\fld2RESULTS\TestDS1_12_12test2\PIC_CL2120 A1 11RUNfinal'
%InputDir='/Users/budoday/Dropbox/WormMovement/AdulanteSwimming/SwimData/fld2RESULTS/TestDS1_12_12test2/'
%InputDir='/Users/budoday/Desktop/SwimData/fld1';

%InputDir='\\murphylab.princeton.edu\data\shared\CL2120_Adult_24h_March_6\with_RNAi';
%'\\murphylab.princeton.edu\data\shared\CL2120_Adult_24h_March_6\with_RNAi'
%pinceton.edu\data\shared\CL2120_Adult_24h_March_6\control_RNAi'
%InputDir='\\murphylab.princeton.edu\data\shared\CL2122_Adult_24h_March_6\with_RNAi_CL2122'
%%InputDir='C:\Documents and Settings\mlab\My Documents\Dropbox\Analysis_Movement'
%Outputdir=InputDir
%trialName='7_4_12Run'

%OneWorm Pipeline
[Alldata]=ProfileWorm(Alldata)
[Alldata]=GetWorm(Alldata)%(InputDir); % 'trialName', trialName
[Alldata]=ProcessSpine(Alldata);%(InputDir);
OneWormFigs(Alldata);

%%TrialStitch one level up!


%MAybe dont use - fix in process spine %%ErrorFixer(InputDir, Outputdir, trialName)
    
%>FancyFIGs
%>OneWormAnalysis

%% Final file checker
%FinalFile =
%'/Users/budoday/Desktop/DATA/Geneva-N2thrash/RESULTS/TestDS1_12_12test2RUNfinal/N2_frame_0179final.mat'
%FinalFile = '/Users/budoday/Desktop/SwimData/fld4RESULTS/TestDS1_12_12test2/PIC_CL2120 A1 5RUNfinal/CL2120 A1 5_frame_0001final'

%FinalFile = '\\murphylab.princeton.edu\data\people\kleemann\WormMovement\AdulanteSwimming\SwimData\fld2RESULTS\TestDS1_12_12test2\PIC_N2_24hRUNfinal\N2 A1 1_frame_0041final'
%OutputChecker(FinalFile)7=======
%OutputChecker(FinalFile)

