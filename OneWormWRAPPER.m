%% WORKFLOW for OneWorm.
% specify your path to the input directory. 
%InputDir='C:\Documents and Settings\mlab\Desktop\reproductive_span_robot\RSbatch2screen1-verify12_2011'
%..InputDir='/Users/budoday/Dropbox/Movement/Movement_Analysis/REDO'

% On a PC always specify a trialName

%%InputDir='C:\Documents and Settings\mlab\My Documents\Dropbox\Analysis_Movement'
%Outputdir=InputDir
%trialName='7_4_12Run'

%%

%OneWorm Pipeline
%[InputDir]=ProfileWorm()
%Outputdir=InputDir;
%[InputDir]=GetWorm(InputDir)%(InputDir); % 'trialName', trialName
%[InputDir]=ProcessSpine(InputDir, Outputdir);%(InputDir);
OneWormFigs(InputDir, Outputdir);
%TrialStitch %stitches together individual replicate (PIC_ folders) and sumarizes data by strain 


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

