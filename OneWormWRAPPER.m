%% WORKFLOW for OneWorm. 
% specify your path in input directory. On a PC
% always specify a trialName


%InputDir='C:\Documents and Settings\mlab\Desktop\reproductive_span_robot\10_28_11_Cheng_Repx1Mated',
%InputDir='C:\Documents and Settings\mlab\Desktop\reproductive_span_robot\RS-Batch1screen1Verif_Prs';
%InputDir='E:\22-Dec-2011-RSrapidVerif'
%InputDir='C:\Documents and Settings\mlab\Desktop\reproductive_span_robot\RSbatch2screen1-verify12_2011'

%InputDir='/Users/budoday/Desktop/DATA/AldulanteSwimmingTest'
%InputDir='/Users/budoday/Desktop/DATA/AldulanteSwimmingTest/N2_24hrun/';
%InputDir='smb://murphylab.princeton.edu/data/people/kleemann/Murphylabdata/Geneva-N2thrash/PICS_A'
%InputDir='/Users/budoday/Desktop/DATA/Geneva-N2thrash'

%InputDir='C:\Documents and Settings\mlab\My Documents\Dropbox\WormMovement\AdulanteSwimming\SwimData\fld2'

%Documents\Dropbox\WormMovement\AdulanteSwimming\SwimData\fld2RESULTS\TestDS1_12_12test2\PIC_CL2120 A1 11RUNfinal'
%InputDir='/Users/budoday/Dropbox/WormMovement/AdulanteSwimming/SwimData/fld2RESULTS/TestDS1_12_12test2/'
%InputDir='/Users/budoday/Desktop/SwimData/fld1';
InputDir='\\murphylab.princeton.edu\data\people\kleemann\WormMovement\AdulanteSwimming\SwimData\fld3'
%InputDir='\\murphylab.princeton.edu\data\people\kleemann\WormMovement\AdulanteSwimming\SwimData\fld1RESULTS\TestDS1_12_12test2'
Outputdir=InputDir;
trialName='TestDS1_12_12test2';

%GetWorm(InputDir, Outputdir, trialName);

InputDir='\\murphylab.princeton.edu\data\people\kleemann\WormMovement\AdulanteSwimming\SwimData\fld2RESULTS\TestDS1_12_12test2';
ProcessSpine(InputDir, Outputdir, trialName);
OneWormFigs(InputDir, Outputdir, trialName)

%% Final file checker

%FinalFile =
%'/Users/budoday/Desktop/DATA/Geneva-N2thrash/RESULTS/TestDS1_12_12test2RUNfinal/N2_frame_0179final.mat'
%FinalFile = '/Users/budoday/Desktop/SwimData/fld4RESULTS/TestDS1_12_12test2/PIC_CL2120 A1 5RUNfinal/CL2120 A1 5_frame_0001final'
%OutputChecker(FinalFile)