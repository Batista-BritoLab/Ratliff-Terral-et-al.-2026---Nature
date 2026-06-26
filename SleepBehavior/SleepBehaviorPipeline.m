%% Sleep behavior pipeline Dreadd Exc - Ratliff, Terral et al., 2026 - Nature


%1/Rename file with session name and extract LFP from Dat
bz_LFPfromDat

%2/ Sleep Scoring
%indicate bad (noisy) channels using 
bz_getSessionInfo(pwd,'editGUI',true);% before running SleepScoreMaster
SleepScoreMaster(pwd,'Notch60Hz',1); 

%2/ Determine duration of time in each state
idx=SleepInStates_ExcDREADD(Pathlist,group);%idx gives sessions indices of Pathlist with restriction, if any
%3/ Display as proportion
StateDistribution_Dreadd(PathlistVeh,PathlistCNO);

%4/ Sleep bouts distribution + NREM onsets - get idx from IndividualMouse_LFP_SOM
SleepBouts_Dreadd(PathlistVeh,PathlistCNO,TitleFig) % TitleFig='Sleep bouts';

%5/ extract low vs high arousal states during wake (QW vs motion)
[LowWAKE,HighWAKE]=GetIntervalLowHighArousalRatliff_Terral2026(basepath);

%6/ extract power spectra in specific states and brain regions
PowerSpec_DeltaAmp(STATE,BrainRegion) % STATE includes WAKE, NREM, REM, LowWAKE, HighWAKE
Figure_LFPDelta_VehvsCNO(PathlistVeh,PathlistCNO,BrainRegion,STATE)%plot powerspec

%7/ Position of the mice in the cage with deeplabcut data
Select_Save_Camera_Channel(Pathlist) %This save the camera channel and create a deeplabcut.mat
MousePosition_display_into_Cage(basepath,Theoframerate)
FiguresDLC %plot figure of proportion, time and distance per state
