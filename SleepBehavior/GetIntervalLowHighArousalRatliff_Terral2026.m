
function [LowWAKE,HighWAKE]=GetIntervalLowHighArousalRatliff_Terral2026(basepath)
%Geoffrey 2026 - see also Ratliff_Terral et al 2026
% Update 7/22/25: change saving path (not anymore FreelingMovingLFP) + load
%the three Acc channels to extract with one is out of range of some time

%To get those intervals, this function used accelerometer (acc) data to
%extract low and high wake activity; save in SleepState.states.mat

sav=1; 
MaxTimeSession=7200;%in sec - usually 2h 
force=1;%if it has already run - will not run if 0

if nargin<1
basepath = pwd; 
else
    cd(basepath)
end

basename = bz_BasenameFromBasepath(basepath);
%% retrieve States

load([basename,'.SleepState.states.mat'])
WAKE=SleepState.ints.WAKEstate;
if force==0
    try
    LowWAKE=SleepState.ints.LowWAKE;
    catch
        return
    end
end
%Restrict to 2h recording
B=SubtractIntervals([0 MaxTimeSession],WAKE);
WAKE=SubtractIntervals([0 MaxTimeSession],B);WAKE(:,3)=WAKE(:,2)-WAKE(:,1);


%% retrieve Acc data
load([basename,'.sessionInfo.mat'])
AccChs=sessionInfo.badchannels(end-2:end);
mChs = AccChs; % zero‑based indices
eegFS = 1250; % sampling rate
AccLFP=bz_GetLFP(mChs,'restrict',[0 MaxTimeSession]);

% N = numel(AccLFP.timestamps);
% mask = false(N,3);
%Restrict mask to wake state
[WAKEindex,~,~] =InIntervals(AccLFP.timestamps,WAKE(:,1:2));%get indexes for wake state
mask=[WAKEindex,WAKEindex,WAKEindex];

for ch=1:3
    A=[];
[periods,~] = Threshold([AccLFP.timestamps,double(AccLFP.data(:,ch))],'<=',0.5e4,'min',0.01);
if ~isempty(periods) && AccumulateTimeInt(periods)>1
A=ConsolidateIntervals(periods,'epsilon',1);
A(A(:,2)-A(:,1)<0.3,:)=[];%exclude short intervals that could be detected but are actual signal
%Now transform in indices and exclude outside if positive deflection before/after becoming crazy
A(:,1)=round(A(:,1)*eegFS-0.5*eegFS);
A(:,2)=round(A(:,2)*eegFS+0.5*eegFS);
end

    if ~isempty(A)
        if A(end,2)>size(mask,1)%if wake cut in the max session time
            A(end,2)=size(mask,1);
        end
        for k = 1:size(A,1)          
            mask(A(k,1):A(k,2),ch) = false;
        end
    end
end

LFP=double(AccLFP.data);
LFP(mask==0)=NaN;
Z = normalize(LFP, 1, "zscore");%1 to do across each column - zscore

%% Calculate local variance using a 3sec window
motion = movvar(sum(abs(Z),2,'omitnan'), 3*eegFS);


% % % %Detect High motion intervals
% % % [Highint, ~] = Threshold([AccLFP.timestamps, motion], '>', 0.5);
% % % 
% % % Highint=ConsolidateIntervals(Highint(:,1:2),'epsilon',5);
% % % %Restrict Intervales to wake (peak can b before because of the movvar)
% % % LowWAKE=SubtractIntervals(WAKE(:,1:2),Highint);
% % % HighWAKE=SubtractIntervals(Highint,LowWAKE);

% Convert WAKE intervals to sample indices
WAKE_idx = (WAKE * eegFS);  % in samples

% Step 2: Create a logical mask of the same length as Accelerometers
mask = false(size(motion));
for i = 1:size(WAKE_idx, 1)
    mask(WAKE_idx(i, 1):WAKE_idx(i, 2)) = true;
end

% Get the Accelerometers signal only within WAKE
restricted_signal = motion(mask);

% Step 4: Compute 80th percentile from restricted signal
thresh = prctile(restricted_signal, 80);  % adjust as needed

% Step 5: Find timepoints that are BOTH in WAKE and above threshold
final_mask = mask & (motion > thresh);

% Step 6: Find start and end indices of high activity
d = diff([0; final_mask(:); 0]);
starts = find(d == 1);
ends   = find(d == -1) - 1;

% Step 7: Convert to time
start_times = (starts - 1) / eegFS;
end_times   = (ends - 1) / eegFS;

% Step 8: Remove short intervals
min_duration = 1;  % in seconds
durations = end_times - start_times;
valid = durations >= min_duration;

Highint = [start_times(valid), end_times(valid)];

% Step 9: Consolidate intervals with gaps < 5 seconds
Highint = ConsolidateIntervals(Highint, 'epsilon', 5);

LowWAKE=SubtractIntervals(WAKE(:,1:2),Highint);
HighWAKE=SubtractIntervals(Highint,LowWAKE);


figure
plot(AccLFP.timestamps,motion,'r')
hold on
PlotIntervals(HighWAKE(:,1:2),'color','k')
hold on
PlotIntervals(LowWAKE(:,1:2),'color','b')
ylabel('Motion (Accel.) LowWAKE (b) and HighWAKE (k)')
xlabel('Time (sec)')
title(basename)
box off

%% Saving
SleepState.ints.LowWAKE=LowWAKE;
SleepState.ints.HighWAKE=HighWAKE;
Accelerometers=motion;

if sav==1
save([basename + ".SleepState.states.mat"],'SleepState')
save([basename + ".AccelerometerSignal.mat"],'Accelerometers')
end

end