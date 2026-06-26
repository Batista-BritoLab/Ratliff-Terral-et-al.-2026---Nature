%Geoffrey 2026 - see also Ratliff, Terral et al 2026

% This script will extract powerspec in a specific state while excluding outlier period. 
% It will also provide delta amplitude by extracting the wavelet 

% At the end of the script, commented, display for each channel difference
% in powerspec between Veh and CNO treatment
%Geoffrey 01/29/24

%Example how to run it:
% PathlistVeh={'E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day2_5610_Saline_230712_093701','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day4_5610_Saline_230714_093205','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day1_5621_Saline_230711_113648','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day3_5621_Saline_230713_112609','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day1_5887_Saline_230711_113611','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day3_5887_Saline_230713_112550','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day2_5890_Saline_230712_135754','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day4_5890_Saline_230714_141716','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day5_5671_Saline_230721_115943','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day7_5671_Saline_230723_113231','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day6_5891_Saline_230722_135441','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day9_5891_Saline_230725_135140','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day6_6211_Saline_230722_092530','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day8_6211_Saline_230724_091819','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day5_6212_Saline_230729_095840','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day7_6212_Saline_230731_100427','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day1_5610_Saline_230924_090941','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day5_5610_Saline_230928_132325','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\hab2_5621_Saline_230923_113204','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\Day4_5621_Saline_230927_112107','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day1_5671_Saline_230924_132048','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day3_5671_Saline_230926_133847','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day2_5887_Saline_230925_111746','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day4_5887_Saline_230927_112114','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day1_6211_Saline_230924_132123','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day3_6211_Saline_230926_133947','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day1_6212_Saline_230924_090625','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day3_6212_Saline_230926_090028'};
% PathlistCNO={'E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day1_5610_CNO_05mgkg_230711_092632','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day3_5610_CNO_05mgkg_230713_091601','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day2_5621_CNO_05mgkg_230712_114729','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day4_5621_CNO_05mgkg_230714_115242','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day2_5887_CNO_05mgkg_230712_114709','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day4_5887_CNO_05mgkg_230714_115227','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day1_5890_CNO_05mgkg_230711_135525','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day3_5890_CNO_05mgkg_230713_133802','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day6_5671_CNO_05mgkg_230722_114550','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day8_5671_CNO_05mgkg_230724_113618','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day5_5891_CNO_05mgkg_230721_141157','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day8_5891_CNO_05mgkg_230724_134757','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day5_6211_CNO_05mgkg_230721_093448','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day7_6211_CNO_05mgkg_230723_091932','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day6_6212_CNO_05mgkg_230730_094227','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day8_6212_CNO_05mgkg_230801_100938','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day2_5610_CNO_05mgkg_230925_090543_unpluglastMin','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day4_5610_CNO_05mgkg_230927_091053','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\Day1_5621_CNO_05mgkg_230924_111709','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\Day3_5621_CNO_05mgkg_230926_113004','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day2_5671_CNO_05mgkg_230925_133729','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day4_5671_CNO_05mgkg_230927_133507','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day1_5887_CNO_05mgkg_230924_111451','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day3_5887_CNO_05mgkg_230926_112730','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day2_6211_CNO_05mgkg_230925_133810','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day4_6211_CNO_05mgkg_230927_133543','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day2_6212_CNO_05mgkg_230925_090228','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day4_6212_CNO_05mgkg_230927_090703'};
% for s=1: numel(PathlistVeh)
%     cd(PathlistVeh{s})
%     PowerSpec_DeltaAmp('NREM','Cortex')
% end

function PowerSpec_DeltaAmp(STATE,BrainRegion)

% STATE='NREM';
MaxTime=120;%in min - max time of the session or analyzis
minTPSRestrict=10;%how long min intervals in sec
THfactor=10;% X times STD - outlier threshold for delta amplitude (to exclude intervals)
passband=[1 4]; %Delta - same as Jacob changed 8/6/2025 from 0.5
% BrainRegion='Cortex';%or 'HPC'
% BrainRegion='HPC';
scalingFactor = 0.195; % Microvolts per raw ADC unit

Originalbasepath = pwd; basename = bz_BasenameFromBasepath(Originalbasepath);

load([basename,'.SleepState.states.mat'])
load([basename,'.sessionInfo.mat'])

% try
%    Channel=sessionInfo.([sprintf(BrainRegion),'Channel']); 
% catch

cd ..
basepath= pwd; base = bz_BasenameFromBasepath(basepath);
load([base,'.Parameters.mat'])
Channel=Parameters.([sprintf(BrainRegion),'Channel']);
cd(Originalbasepath)
% end


if strcmp(STATE,'All')
    Interval=[0 SleepState.idx.timestamps(end,1)];
else
StateIntervals=SleepState.ints;
try
Interval=StateIntervals.([sprintf(STATE),'state']);
catch
Interval=StateIntervals.([sprintf(STATE)]);   
end
end

    if strcmp(STATE,'LowWAKE') || strcmp(STATE,'HighWAKE')
        A=SubtractIntervals([20*60 60*MaxTime],Interval);%exclusion of first 20min to let drug to act
        Interval=SubtractIntervals([20*60 60*MaxTime],A);

    else

        A=SubtractIntervals([0 60*MaxTime],Interval);
        Interval=SubtractIntervals([0 60*MaxTime],A);
    end

try
Interval(:,3)=Interval(:,2)-Interval(:,1);%duration
Interval_min_tps=Interval;
Interval_min_tps(Interval_min_tps(:,3)<minTPSRestrict,:)=[]; %do not consider if interval shorter than Xsec
catch
 Interval_min_tps=[];
end

if isempty(Interval_min_tps)|| AccumulateTimeInt(Interval_min_tps)<=40
    alreadysaved=true;
    eval(sprintf('%s.%s = "[]";', STATE,sprintf('Power')))
else
[lfp] = bz_GetLFP(sessionInfo.AnatGrps(1).Channels,'restrict',Interval_min_tps); %in sec

%Extract LFP
x1=[];lfpTps=[];
for i=1:size(Interval_min_tps,1)
x=double(lfp(i).data);
x1=[x1; x];
lfpTps=[lfpTps;(lfp(i).timestamps)];
end
Median=median(x1,2); %median references over time
C=find(lfp(1).channels==Channel);
MedianRef=x1(:,C)-Median;
MedianRef=MedianRef*scalingFactor;

%Calcul Delta Power Amplitude
[wave,f,t,coh,wphases,raw,coi,scale, period, scalef]=getWavelet(MedianRef,1250,passband(1),passband(2),8,0); % 


%Exclude outliers based on amplitude signal
% A=max(wave);%change 09/18/23: we don't need to specify what freq, we want all passband range
A=sum(wave);
STD1 = std(A);
TH = THfactor*STD1;
[periods,in1] = Threshold([lfpTps,A'],'>',TH);
            if max(in1)
                AClean = A(~in1);
                MedianRef = MedianRef(~in1);
                lfpTps = lfpTps(~in1);

                LFPTPS=nan(size(in1,1),1);%ADD 09/06/23
                LFPTPS(in1~=1)=lfpTps;
                lfpTps=LFPTPS;
                
%               figure
%               plot(lfpTpsClean,AClean)
    
                Interval_min_tps=SubtractIntervals(Interval_min_tps(:,1:2),periods);

                Interval_min_tps(:,3)=Interval_min_tps(:,2)-Interval_min_tps(:,1);
%                 Interval_min_tps(Interval_min_tps(:,3)<minTPSRestrict,:)=[]; %do not consider if interval shorter than 10sec

                [Power, freq]=pwelch(MedianRef,[1250],[],[12500],1250,'onesided','power');
                [wave,f,t,coh,wphases,raw,coi,scale, period, scalef]=getWavelet(MedianRef,1250,passband(1),passband(2),8,0); % 
%                 A=max(wave);
                A=sum(wave);
                Amp=nan(size(in1,1),1);%ADD 09/06/23
                Amp(in1~=1)=A;
                A=Amp;
            else
                [Power, freq]=pwelch(MedianRef,[1250],[],[12500],1250,'onesided','power');
            end

% Create the dynamic field name
AmplitudeDelta=nanmean(A);
Power(1)=[];
freq(1)=[];

%saving
try
    load([basename,'.powerLFP.mat'])
    alreadysaved=1;
catch
    alreadysaved=false;
end

eval(sprintf('%s.%s = Power;', STATE,sprintf('%sPowerSpecData', BrainRegion)))
eval(sprintf('%s.%s = freq;', STATE,sprintf('%sPowerSpecFreq', BrainRegion)))
eval(sprintf('%s.%s = "Pwelch";', STATE,sprintf('Method')))
eval(sprintf('%s.%s = Interval_min_tps;', STATE,sprintf('%sInterval_min_tps', BrainRegion)))%before it was Interval instead of Interval_min_tps - 07/10/2025
eval(sprintf('%s.%s = Channel;', STATE,sprintf('%sChannel', BrainRegion)))
eval(sprintf('%s.%s = Channel;', STATE,sprintf('%sChannel', BrainRegion)))
eval(sprintf('%s.%s = MedianRef;', STATE,sprintf('%sLFPMedianRef', BrainRegion)))
eval(sprintf('%s.%s = AmplitudeDelta;', STATE,sprintf('%sMeanAmplitudeDelta', BrainRegion)))


% PowerSpec.(sprintf(STATE)).([sprintf(BrainRegion),'PowerSpecData'])=Power;
% PowerSpec.(sprintf(STATE)).([sprintf(BrainRegion),'PowerSpecFreq'])=freq;
% PowerSpec.(sprintf(STATE)).Method='Pwelch';
% PowerSpec.(sprintf(STATE)).([sprintf(BrainRegion),'Interval'])=Interval_min_tps;
% PowerSpec.(sprintf(STATE)).([sprintf(BrainRegion),'Channel'])=Channel;
% % DeltaAmp.(sprintf(STATE)).([sprintf(BrainRegion),
end

if alreadysaved
save (basename + ".powerLFP.mat" , STATE,'-v7.3','-append')
else
save (basename + ".powerLFP.mat" , STATE,'-v7.3')
end

end


%% 
% % % %% DISPLAY difference of powerspec for all session pairs (CNO-Veh) for all channels
% % % BrainRegion='Cortex';
% % % minTPSRestrict=10;
% % % scalingFactor = 0.195; % Microvolts per raw ADC unit
% % % 
% % % 
% % % PathlistVeh={'E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day2_5610_Saline_230712_093701','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day4_5610_Saline_230714_093205','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day1_5621_Saline_230711_113648','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day3_5621_Saline_230713_112609','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day1_5887_Saline_230711_113611','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day3_5887_Saline_230713_112550','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day2_5890_Saline_230712_135754','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day4_5890_Saline_230714_141716','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day5_5671_Saline_230721_115943','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day7_5671_Saline_230723_113231','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day6_5891_Saline_230722_135441','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day9_5891_Saline_230725_135140','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day6_6211_Saline_230722_092530','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day8_6211_Saline_230724_091819','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day5_6212_Saline_230729_095840','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day7_6212_Saline_230731_100427','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day1_5610_Saline_230924_090941','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day5_5610_Saline_230928_132325','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\hab2_5621_Saline_230923_113204','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\Day4_5621_Saline_230927_112107','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day1_5671_Saline_230924_132048','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day3_5671_Saline_230926_133847','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day2_5887_Saline_230925_111746','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day4_5887_Saline_230927_112114','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day1_6211_Saline_230924_132123','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day3_6211_Saline_230926_133947','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day1_6212_Saline_230924_090625','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day3_6212_Saline_230926_090028'};
% % % PathlistCNO={'E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day1_5610_CNO_05mgkg_230711_092632','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Day3_5610_CNO_05mgkg_230713_091601','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day2_5621_CNO_05mgkg_230712_114729','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Day4_5621_CNO_05mgkg_230714_115242','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day2_5887_CNO_05mgkg_230712_114709','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Day4_5887_CNO_05mgkg_230714_115227','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day1_5890_CNO_05mgkg_230711_135525','E:\Ephy\homeCage\CNO_0.5mgkg\5890\Day3_5890_CNO_05mgkg_230713_133802','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day6_5671_CNO_05mgkg_230722_114550','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Day8_5671_CNO_05mgkg_230724_113618','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day5_5891_CNO_05mgkg_230721_141157','E:\Ephy\homeCage\CNO_0.5mgkg\5891\Day8_5891_CNO_05mgkg_230724_134757','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day5_6211_CNO_05mgkg_230721_093448','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Day7_6211_CNO_05mgkg_230723_091932','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day6_6212_CNO_05mgkg_230730_094227','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Day8_6212_CNO_05mgkg_230801_100938','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day2_5610_CNO_05mgkg_230925_090543_unpluglastMin','E:\Ephy\homeCage\CNO_0.5mgkg\5610\Round2\Day4_5610_CNO_05mgkg_230927_091053','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\Day1_5621_CNO_05mgkg_230924_111709','E:\Ephy\homeCage\CNO_0.5mgkg\5621\Round2\Day3_5621_CNO_05mgkg_230926_113004','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day2_5671_CNO_05mgkg_230925_133729','E:\Ephy\homeCage\CNO_0.5mgkg\5671\Round2\Day4_5671_CNO_05mgkg_230927_133507','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day1_5887_CNO_05mgkg_230924_111451','E:\Ephy\homeCage\CNO_0.5mgkg\5887\Round2\Day3_5887_CNO_05mgkg_230926_112730','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day2_6211_CNO_05mgkg_230925_133810','E:\Ephy\homeCage\CNO_0.5mgkg\6211\Round2\Day4_6211_CNO_05mgkg_230927_133543','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day2_6212_CNO_05mgkg_230925_090228','E:\Ephy\homeCage\CNO_0.5mgkg\6212\Round2\Day4_6212_CNO_05mgkg_230927_090703'};
% % % 
% % % PairSessions=[ [1:2:numel(PathlistVeh)]' [2:2:numel(PathlistVeh)]'];
% % % 
% % % for K=1:size(PairSessions,1)
% % %     PathlistV=PathlistVeh(PairSessions(K,:));
% % %     PathlistC=PathlistCNO(PairSessions(K,:));
% % % 
% % % figure
% % % for s=1:numel(PathlistV)
% % %     k=s;
% % %     %retrieve Data for Veh
% % % cd(PathlistV{s})
% % % basepath = pwd; basename = bz_BasenameFromBasepath(basepath);
% % % load(basename + ".FreelyMovingLFP" + minTPSRestrict+ sprintf(BrainRegion)+ ".mat")
% % % load([basename,'.sessionInfo.mat'])
% % % ChannelsVeh = sessionInfo.AnatGrps(1).Channels;%exclude bad channels
% % % Interval_min_tpsVeh=FreelyMovingLFP.(sprintf(STATE)).Power.Interval(:,1:2);
% % % [lfp1] = bz_GetLFP(sessionInfo.AnatGrps(1).Channels,'restrict',Interval_min_tpsVeh); %in sec
% % % 
% % % %Extract LFP
% % % x1=[];
% % % for i=1:size(Interval_min_tpsVeh,1)
% % % x=double(lfp1(i).data);
% % % x1=[x1; x];
% % % end
% % % MedianVeh=median(x1,2); %median references over time
% % % 
% % %     %retrieve Data for CNO
% % % cd(PathlistC{s})
% % % basepath = pwd; basename = bz_BasenameFromBasepath(basepath);
% % % load(basename + ".FreelyMovingLFP" + minTPSRestrict+ sprintf(BrainRegion)+ ".mat")
% % % load([basename,'.sessionInfo.mat'])
% % % ChannelsCNO = sessionInfo.AnatGrps(1).Channels;%exclude bad channels
% % % 
% % % Interval_min_tpsCNO=FreelyMovingLFP.(sprintf(STATE)).Power.Interval(:,1:2);
% % % [lfp2] = bz_GetLFP(sessionInfo.AnatGrps(1).Channels,'restrict',Interval_min_tpsCNO); %in sec
% % % 
% % % %Extract LFP
% % % x2=[];
% % % for i=1:size(Interval_min_tpsCNO,1)
% % % x=double(lfp2(i).data);
% % % x2=[x2; x];
% % % end
% % % MedianCNO=median(x2,2); %median references over time
% % % 
% % % Channels=ChannelsVeh(ChannelsVeh==ChannelsCNO);%keep only common channel
% % % NCh=numel(Channels);
% % % 
% % % cd ..
% % % basepath = pwd; basename = bz_BasenameFromBasepath(basepath);
% % % load([basename,'.Parameters.mat'])
% % % CortexChannel=Parameters.CortexChannel;
% % % 
% % % for nCh=1:NCh
% % % MedianRefVeh=[];
% % % C=find(lfp1(1).channels==Channels(nCh));
% % % MedianRefVeh=x1(:,C)-MedianVeh;
% % % C=find(lfp2(1).channels==Channels(nCh));
% % % MedianRefCNO=x2(:,C)-MedianCNO;
% % % 
% % % %Calcul Delta Power
% % % [PowerVeh(:,nCh), freq]=pwelch(MedianRefVeh*scalingFactor,[1250],[],[],1250,'onesided','power');
% % % [PowerCNO(:,nCh), freq]=pwelch(MedianRefCNO*scalingFactor,[1250],[],[],1250,'onesided','power');
% % % 
% % % subplot(NCh,numel(PathlistV),k)
% % % if Channels(nCh)==CortexChannel
% % %     COLOR='r';
% % % else
% % %     COLOR='b';
% % % end
% % % plot(freq(1:50),PowerCNO(1:50,nCh)-PowerVeh(1:50,nCh),COLOR)
% % % box off
% % % title(['Channel ',num2str(Channels(nCh))])
% % % k=k+2;
% % % if nCh==NCh
% % %     xlabel('Freq. (Hz)')
% % % end
% % % end
% % % sgtitle(basename)
% % % end
% % % 
% % % end
