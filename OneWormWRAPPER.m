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
InputDir='/Users/budoday/Desktop/DATA/Geneva-N2thrash'
Outputdir=InputDir;
trialName='TestDS1_12_12test2';

%GetWorm(InputDir, Outputdir, trialName)
OneWormFigs(InputDir, Outputdir, trialName)