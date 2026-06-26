%Geoffrey 2026 - see Ratliff_Terral et al 2026

function idx=SleepInStates_ExcDREADD(Pathlist,group)
%Change in % time instead of raw min
% Pathlist=[PathlistVeh,PathlistCNO];
Data_per_mouse=1;%1 to generate data per mouse, else it will be per session. If per mouse, it assumes that 2 sessions/mouse successively presented in Pathlist
% Data_per_mouse=0;
DATA=[];%to get all the data into one variable

Veh_WAKE='short';%restrict analysis to short veh sessions >60% WAKE out of 2h - paired veh-cno restriction
% Veh_WAKE='long';% else = no restriction - all 2h session
% Veh_WAKE=0;
%% % % % group{1}=[ones(8,1)*1;ones(8,1)*2;ones(8,1)*3]; %%How many mice
if nargin<2
 group{1}=[ones(numel(Pathlist),1)*1];
end


% Within each state interval; column 1-2 start stop interval, column 3
% duration, column 4 cumulative time. So last value of 4th colum is max time in this state

force=0;
MaxTime=120;%min

STATEintervalWNR={'WAKE','NREM','REM'};
% STATEintervalWNR={'WAKE','LowWAKE','HighWAKE'}; %restrict only after 20min - same as for power
% 
% group{1}=[2;1;2;1;2;1;2;1];%group 1 is the treatment: 1 veh vs 2 drug
% group{1}=[1;2;1;2;1;2;1;2];%group 1 is the treatment: 1 veh vs 2 drug
% 
% group{2}=1;%group 2 is the mouse;1 if only one mouse; 2 if 2 mice...
% group{2}=ones(numel(group{1}),1)*group{2};
% % group {2}=ones(8,1)*1;
% % group{2}=[1;1;1;1;2;2;2;2];
% 
COLORCODE={'k','r'};

%% Find indices according to exp condition (veh vs drug vs what mouse)

TREATMENT=[];
for s=1:numel(Pathlist)   
    cd(Pathlist{s})
    basepath = pwd; basename = bz_BasenameFromBasepath(basepath);
    load([basename,'.SleepState.states.mat'])
    StateIntervals=SleepState.ints;
    load([basename,'.sessionInfo.mat'])
    Treatment=sessionInfo.Treatment;
    
    for i=1:3 %3 states: WAKE, NREM, REM
        Interval=[];
    STATE=STATEintervalWNR{1,i};
    try
    Interval=StateIntervals.([sprintf(STATE),'state']);
    catch
        Interval=StateIntervals.([sprintf(STATE)]);
        % load([basename,'.FreelyMovingLFPCortex.mat'])
        % Interval=FreelyMovingLFP.(sprintf(STATE)).Interval(:,1:2);
    end
       
    if strcmp(STATEintervalWNR{1,2}, 'LowWAKE') && i>1
    A=SubtractIntervals([20*60 60*MaxTime],Interval);%restrict only after 20min - same as for power
    Interval=SubtractIntervals([20*60 60*MaxTime],A);
    else
    A=SubtractIntervals([0 60*MaxTime],Interval);
    Interval=SubtractIntervals([0 60*MaxTime],A);
    end
    Interval(:,3)=Interval(:,2)-Interval(:,1);%duration
    Interval(:,4)=cumsum(Interval(:,3));
    STATEintervalWNR{1+s,i}=Interval;
    end
    if strcmp(Treatment,'Vehicle')
        TREATMENT=[TREATMENT,1];
    elseif strcmp(Treatment,'J60')|| strcmp(Treatment,'CNO')
        TREATMENT=[TREATMENT,2];
    end
end
  
DrugMouseIDX=[];VehMouseIDX=[];
for ii=1:max(group{1})
    for jj=1:max(TREATMENT)
        if jj==1 
           VehMouseIDX{1,ii}=find(TREATMENT'==jj & group{1}==ii);
        elseif jj==2
           DrugMouseIDX{1,ii}=find(TREATMENT'==jj & group{1}==ii);
        end
    end
end


fig=figure;
kk=0;
for k=1:max(group{1})%removing 2 groups - here only how many mice so group{2} initially is now group{1}
% % for k=1:max(group{2})
    %Raw time
    for i=1:3 %3 states: WAKE, NREM, REM
    STATE=STATEintervalWNR{1,i};

    NotemptyCells = find(~cellfun(@isempty, STATEintervalWNR(2:end, i)));
    emptyCells = find(cellfun(@isempty, STATEintervalWNR(2:end, i)));
    CumTimeInt=[];
    CumTimeInt(emptyCells,1)=[0];
    CumTimeInt(NotemptyCells,1)=(cellfun(@(c) c(end,4), STATEintervalWNR(1+NotemptyCells, i)));
    
    subplot(max(group{1}),3,i+kk)
% %     subplot(max(group{2}),3,i+kk)
    % COLORCODE={'r','k'};%inverted for boxplot
    COLORCODE={'k','r'};
    conditions={'Veh','CNO'};
    ytitle={['Time in ',STATE,' (%)']};
    
    data=[];
    data{1}=CumTimeInt(VehMouseIDX{1,k})./60;
    data{2}=CumTimeInt(DrugMouseIDX{1,k})./60;

    if Data_per_mouse ==1
    MMM=[];m=0;JJJ=[];
    Nmice=numel(VehMouseIDX{1})/2;%7/7/2025
        for M=1:Nmice
            MM=[];JJ=[];
            MM=(data{1}(1+m)+data{1}(2+m))/2;
            JJ=(data{2}(1+m)+data{2}(2+m))/2;
            MMM=[MMM;MM];JJJ=[JJJ;JJ];
            m=m+2;
        end
        data{1}=MMM;
        data{2}=JJJ;
    end
    %restrict analysis only with veh sessions short or long
    if strcmp(Veh_WAKE,'short')
        if i==1
        idx=(data{1}>48);%120-120*0.6
        disp(['Restriction N: ',num2str(sum(idx)),' out of ', num2str(size(data{1},1))]);
        end
    elseif strcmp(Veh_WAKE,'long')
        if i==1
        idx=(data{1}<=48);
        disp(['Restriction N: ',num2str(sum(idx)),' out of ', num2str(size(data{1},1))]);
        end
    else
        idx=logical(ones(size(data{1})));
    end
    data{1}=data{1}(idx)./MaxTime*100; %change min to %
    data{2}=data{2}(idx)./MaxTime*100;

    %BAR GRAPH - it works if not same size of data
    [~,datamean,dataSEM,pval]=MeanSEMgraph(data, conditions, ytitle, COLORCODE,'paired');

    s1 = std(data{1});
    s2 = std(data{2});
    n1 = length(data{1});
    n2 = length(data{2});

    pooled_std = sqrt(((n1 - 1)*s1^2 + (n2 - 1)*s2^2) / (n1 + n2 - 2));
    d = (mean(data{2})-mean(data{1})) / pooled_std;%calculate size effect
    disp(['size_effect ',num2str(d)])

    if i==3%REM
        if strcmp(STATEintervalWNR{1,2},'LowWAKE')
            ylim([0 60])
        else
            ylim([0 10])
        end
    else
    ylim([0 100])
    end

    hold on
    DATA=[DATA,data{1},data{2}];
    end
kk=kk+max(group{1});
% % kk=kk+max(group{2});
end
%     %% Cumulative
%     AllIntervals=0:IntDuration:MaxTime*60;

set(fig, 'Position', [500, 500, 900, 500]);
if strcmp(STATEintervalWNR{1,2},'LowWAKE')
figure
% [~,datamean,dataSEM,pval]=MeanSEMgraph([DATA(:,3)./DATA(:,1),DATA(:,4)./DATA(:,2)], conditions, 'Proportion Low WAKE / total WAKE', COLORCODE,'paired');
[~,datamean,dataSEM,pval]=MeanSEMgraph([DATA(:,5)./DATA(:,3),DATA(:,6)./DATA(:,4)], conditions, 'Proportion High WAKE / Low WAKE', COLORCODE,'paired');
ylim([0 0.4])
end
end