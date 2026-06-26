function stateInts = getStatesIntervals(basepath,facemotion, facemotionthreshold,IntanOnset)
%dependency: defMovementIntervals and extractArousalState, and of course to have sleep scoring and facemotion

%Careful: facemotion time was already aligned to imaging, add this time(1)
%to get same alignement as intan!
basename = bz_BasenameFromBasepath(basepath);

try
load(fullfile(basepath,basename + ".IntStates.mat"),"stateInts")
catch

% Load facemotion information
%% facemotion = extractArousalState(basepath);%Jacob
% %Alternative to not rerun extractArousal but still taking raw motion_1
% facemapOutputFile = dir(fullfile(string(basepath),'*_proc.mat'));%Geoffrey
% load(fullfile(facemapOutputFile.folder, facemapOutputFile.name), 'motion_1');
% 
% fsfacemotion=1/mean(diff(facemotion.time));
% smoothingKernel = 4*fsfacemotion; % samples
% facemotion.motionSmooth = smoothdata(motion_1, 'gaussian', smoothingKernel);%use original signal - same as Jacob
% ntrigs = length(facemotion.time);
% queriedValues = rescale(1:ntrigs, 1, length(motion_1));
% facialMotion = interp1(1:length(motion_1),motion_1, queriedValues);
% facemotion.motionSmooth = smoothdata(facialMotion, 'gaussian', smoothingKernel);
%%

% smoothingKernel = 100; % samples

% % smooth signal
fsfacemotion=1/mean(diff(facemotion.time));
smoothingKernel = 4*fsfacemotion; % samples
motionSmooth = smoothdata(facemotion.signal, 'gaussian', smoothingKernel);

%estimation of baseline and normlaize signal
x_pos = motionSmooth(motionSmooth >= 0);
if isempty(x_pos)
    baseline = median(motionSmooth);
else
    baseline = prctile(x_pos, 5); %the 5th percentile of non-negative values as baseline
end
%subtract baseline
facemotion.motionSmooth = motionSmooth - baseline;

facemotion.time=facemotion.time-IntanOnset;%to re aligned to intan

figure
plot(facemotion.time,facemotion.motionSmooth,'k'); ylim([0 20]); hold on; yline(facemotionthreshold,'r','LineWidth',2); ylabel('Facemotion')

% Get state intervals
load(fullfile(basepath,basename + ".SleepState.states.mat"),"SleepState")
stateInts.SWS = SleepState.ints.NREMstate;
stateInts.REM = SleepState.ints.REMstate;
stateInts.WAKE = SleepState.ints.WAKEstate;
stateInts.Movement = facemotion.time(defMovementIntervals(facemotion.motionSmooth', facemotionthreshold, 1));
% [~, sleepMask] = RestrictInts(stateInts.Movement, sort([stateInts.REM; stateInts.SWS]));
stateInts.Movement=SubtractIntervals(stateInts.Movement, sort([stateInts.REM; stateInts.SWS]));
% stateInts.Movement(sleepMask,:) = [];
stateInts.stillInts = [stateInts.Movement(1:end-1,2) stateInts.Movement(2:end,1)];
stateInts.stillInts = stateInts.stillInts(diff(stateInts.stillInts,1,2) > 5,:); % For stillness
stateInts.Movement = stateInts.Movement(diff(stateInts.Movement,1,2) > 5,:);
stateInts.QuietWAKE = IntersectIntervals(stateInts.WAKE, stateInts.stillInts); % For stillness w/o sleep



fullFilePath=fullfile(basepath,[basename, '.IntStates.mat']);
save(fullFilePath, 'stateInts')
end
end