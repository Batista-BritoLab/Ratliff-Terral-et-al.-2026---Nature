%% 2p pipeline - Ratliff, Terral et al., 2026 - Nature

%% Retrieve 2p signal and save output signals
s2pOut = retrieveSuite2pOutput(imgpath + "\suite2p",[]);%if no pre img folder
fs=s2pOut.ops.fs;
ops.normalizationtype='MovMedian';
%create data structure
ImgStruct = initializeImagingAnalysisStruct(s2pOut,fs,ops);
imgStruct = align2pData_Pre_Post(ImgStruct,[basepath+ '\' + basename + '.dat'] ,scopeTriggerChannel, nChans, fsAnalog, treadmillChan, cameraTriggerChan);
save (basename + ".ImgStruct.mat", 'imgStruct', '-v7.3')%save in individual path if several files
%Save output in common path - compile pre and post recording if any
Aligned_Motion_Pre_Post


%% retrict calcium activity to states

facemotionthreshold=1000;
facemotionthreshold=1;
ActivitySTATES=[];

for s=1:numel(Pathlist)%list of pathlist for individual sessions
    cd(Pathlist{s})
    basename=bz_BasenameFromBasepath(Pathlist{s})
    load([basename,'.outputStruct.mat'])
   
    stateInts = getStatesIntervals(Pathlist{s}, facemotion,facemotionthreshold,treadmill.time(1));

    CalciumStates=RestrictCalciumtoStates(img,stateInts,treadmill.time(1));
    [ranova_results,BonfStat]=figureCalciumActivity_inStates(CalciumStates, 'DF/F');%quick plot

ActivitySTATES=[ActivitySTATES;CalciumStates];
end

%% calcium activity correlation to state measurements
State_metrics_vs_cell_activity_correlation_Ratliff_Terral2026

%% calcium activity vs up/down states with flexible probes
Visualize_DOWN_states_Ratliff_Terral2026 %detect UP/DOWN states
Sync_Imaging_down_Ratliff_Terral2026 %synchronize signals
MeanCa_UP_DOWN_UP_Ratliff_Terral2026%plot mean Ca+ according to UP and DOWN states