function SleepBouts_Dreadd(PathlistVeh,PathlistCNO,TitleFig)
%Geoffrey 3/31/2026: adapted from Ratliff*,Terral* 2026
% % SleepBouts_Dreadd(PathlistVeh(repelem(idx,2)),PathlistCNO(repelem(idx,2)),TitleFig)
%Update 7/17/15 to not sue FreelyMovingLFP but direct state in SleepState -
%default use 2h for analysis but can change MaxTime otherwise + add
%analysis per mouse Data_per_mouse=1;%1 to generate data per mouse, else it will be per session. If per mouse, it assumes that 2 sessions/mouse successively presented in Pathlist
%add idx if restriction of mice

% cd("E:\Ephy\homeCage\CNO_0.5mgkg")
% load('Parameters.mat')
% PathlistVeh=Inh_Dreadd.PathlistVeh;
% PathlistCNO=Inh_Dreadd.PathlistCNO;
% PathlistVeh=Exc_Dreadd.PathlistVeh;
% PathlistCNO=Exc_Dreadd.PathlistCNO;
%Here we calculate for each bin of 30sec how many sleeps intervals occur as
%compared to all sleep intervals

%Also calculate Sleep Onsets
MaxTime=120;%in min - max time of the session or analyzis
Data_per_mouse=1;%1 to generate data per mouse, else it will be per session. If per mouse, it assumes that 2 sessions/mouse successively presented in Pathlist

% TitleFig='CtlDreadd';
% TitleFig='InhDreadd';
% TitleFig='ExcDreadd';

v=0;j=0;IntervalVeh=[];IntervalDRUG=[];
STATE='NREM';
BrainRegion='Cortex';
Pathlist=[PathlistVeh,PathlistCNO];
for s=1:numel(Pathlist)
    cd(Pathlist{s})
INT=[];
basename = bz_BasenameFromBasepath(pwd);
% try
% load(basename+".FreelyMovingLFP10"+sprintf(BrainRegion)+".mat")
% catch
% load(basename+".FreelyMovingLFP"+sprintf(BrainRegion)+".mat")    
% end

load([basename + ".SleepState.states.mat"])
load([basename + ".sessionInfo.mat"])

    if strcmp(sessionInfo.Treatment,'Vehicle')
        v=v+1;
        INT=SleepState.ints.([sprintf(STATE),'state']);
        A=SubtractIntervals(INT,[0 MaxTime*60]);%to extract only 2h
        INT=SubtractIntervals(INT,A);
        IntervalVeh{v}=[INT,INT(:,2)-INT(:,1)];
       % IntervalVeh{v}=FreelyMovingLFP.(sprintf(STATE)).Interval;
    elseif strcmp(sessionInfo.Treatment,'J60')||strcmp(sessionInfo.Treatment,'CNO')
        j=j+1;
        INT=SleepState.ints.([sprintf(STATE),'state']);
        A=SubtractIntervals(INT,[0 MaxTime*60]);%to extract only 2h
        INT=SubtractIntervals(INT,A);
        IntervalDRUG{j}=[INT,INT(:,2)-INT(:,1)];
        % IntervalDRUG{j}=FreelyMovingLFP.(sprintf(STATE)).Interval;
    end

end

colorcodeV='k';
colorcodeCNO='r';


%% Distribution of sleep events
IndivSession=0;
% % Individual sessions:
%Time course of sleep bouts hist
if IndivSession==1
figure
for s=1:numel(IntervalVeh)
subplot(1,numel(IntervalVeh),s)

hold on
bar(IntervalVeh{s}(:,3),'EdgeColor', 'none','FaceColor','k','FaceAlpha',0.5)
box off
ylabel('Bout duration (sec)')
xlabel('Chronological bout occurence')
end

figure
for s=1:numel(IntervalDRUG)
subplot(1,numel(IntervalDRUG),s)

hold on
bar(IntervalDRUG{s}(:,3),'EdgeColor', 'none','FaceColor','r','FaceAlpha',0.5)
box off
ylabel('Bout duration (sec)')
xlabel('Chronological bout occurence')
end

%Distribution of sleep bouts hist
% figure
% subplot(2,2,4)
% else
% subplot(1,2,2)
% end
% hold on

figure
for s=1:numel(IntervalVeh)

if IntervalVeh{s}(1,3)~=0
[n,edges,bin] = histcounts(IntervalVeh{s}(:,3),'BinWidth',[30]);
bar(edges(1:end-1),[n],'EdgeColor', 'none','FaceColor','k','FaceAlpha',0.5)
else
    bar(IntervalVeh{s}(:,3),'EdgeColor', 'none','FaceColor','k','FaceAlpha',0.5)
    n=0;
end
end
box off
figure
for s=1:numel(IntervalDRUG)

if IntervalDRUG{s}(1,3)~=0
[n,edges,bin] = histcounts(IntervalDRUG{s}(:,3),'BinWidth',[30]);
bar(edges(1:end-1),[n],'EdgeColor', 'none','FaceColor','r','FaceAlpha',0.5)
else
    bar(IntervalDRUG{s}(:,3),'EdgeColor', 'none','FaceColor','r','FaceAlpha',0.5)
    n=0;
end
end
box off
end


%%All sesssions together
%Here we calculate for each bin of 30sec how many sleeps intervals occur as
%compared to all sleep intervals
n=[];
nCNO=[];
% EDGES=0:30:300;
EDGES=0:180:540;%used for initial submission
% EDGES=[0 180 1200];
% EDGES=0:120:480;
EDGES(end)=3600;%I don't think there is an interval above 60min!
Counts=[];

for s=1:numel(IntervalVeh)
[Counts,~,~] = histcounts(IntervalVeh{s}(:,3),EDGES);
n(s,:)=Counts./sum(Counts);%normalized to occurence
if isnan(n(s))%if no sleep then goes to first bout duration
    n(s,:)=[1 zeros(1,size(EDGES,2)-2)];
end
end
Counts=[];
for s=1:numel(IntervalDRUG)
[Counts,~,~] = histcounts(IntervalDRUG{s}(:,3),EDGES);
nCNO(s,:)=Counts./sum(Counts);%normalized to occurence
if isnan(nCNO(s))%if no sleep then goes to first bout duration
    nCNO(s,:)=[1 0 0];
end
end
% figure
% for s=1:numel(IntervalVeh)
% % A=ksdensity(,n(s,:),1:size(EDGES(2:end),2));
% hold on
% plot(n(s,:));
% end



%% to compile data per mouse
    if Data_per_mouse ==1
    MMM=[];m=0;JJJ=[];
    Nmice=size(n,1)/2;
        for M=1:Nmice
            MM=[];JJ=[];
            MM=nanmean(n(1+m:2+m,:),1);
            JJ=nanmean(nCNO(1+m:2+m,:),1);
            MMM=[MMM;MM];JJJ=[JJJ;JJ];
            m=m+2;
        end
        n=MMM;
        nCNO=JJJ;
    end





MEANcounts=mean(n);
SEMcounts=std(n,[],1)./sqrt(size(n,1)); 
MEANcountsCNO=mean(nCNO);
SEMcountsCNO=std(nCNO,[],1)./sqrt(size(nCNO,1)); 

figure
hb1=shadedErrorBar(1:size(EDGES(2:end),2),MEANcounts,SEMcounts);
hb1.mainLine.Color=colorcodeV;
hb1.patch.FaceColor=colorcodeV;
hb1.patch.FaceAlpha=0.5;
set(hb1.edge,'LineWidth',2,'LineStyle','none');
hold on
hb2=shadedErrorBar(1:size(EDGES(2:end),2),MEANcountsCNO,SEMcountsCNO);
hb2.mainLine.Color=colorcodeCNO;
hb2.patch.FaceColor=colorcodeCNO;
hb2.patch.FaceAlpha=0.5;
set(hb2.edge,'LineWidth',2,'LineStyle','none');

% ax=gca;ax.XTick=[1:10];ax.XTickLabel=num2cell(EDGES(2:end));ax.XTickLabel{end}=['>',num2str(EDGES(end-1)+EDGES(2))];
ax=gca;ax.XTick=[1:3];ax.XTickLabel=num2cell(EDGES(2:end));ax.XTickLabel{end}=['>',num2str(EDGES(end-1)+EDGES(2))];

ylabel('Probability')
xlabel('Bouts duration (sec)')
% xlim([0 11])
xlim([0 4])
titleFIG=['SleepBouts_',TitleFig];
%alternative figure with errorbar
figure
errorbar(EDGES(1:end-1),MEANcounts,SEMcounts,colorcodeV,'LineWidth',2)
hold on
errorbar(EDGES(1:end-1),MEANcountsCNO,SEMcountsCNO,colorcodeCNO,'LineWidth',2)
box off
xlim([(EDGES(1))-10 (EDGES(end-1))+10])
ylabel('Probability')
xlabel('Bouts duration (sec)')
% xlim([0 11])
% xlim([0 4])
titleFIG=['SleepBouts_',TitleFig];
title(titleFIG)
% cd('E:\Ephy\homeCage\CNO_0.5mgkg')
% savefig(titleFIG)
ax=gca;ax.XTick=[0,180,360];ax.XTickLabel=num2cell(EDGES(2:end));ax.XTickLabel{end}=['>',num2str(EDGES(end-1)+EDGES(2))];

%Two way anova
T = table;
T.Subject = (1:size(n,1))';            % subject IDs
T.n1 = n(:,1); T.n2 = n(:,2); T.n3 = n(:,3);
T.nC1 = nCNO(:,1); T.nC2 = nCNO(:,2); T.nC3 = nCNO(:,3);

% Convert to long format for anovan
Treatment = [repmat({'n'},3,1); repmat({'nCNO'},3,1)];
TimeBin   = repmat((1:3)',2,1);
within = table(Treatment, TimeBin);
rm = fitrm(T, 'n1-nC3 ~ 1', 'WithinDesign', within);
ranovatbl = ranova(rm, 'WithinModel', 'Treatment*TimeBin');

% Bonferroni correction
pvals = nan(1,3);
for c = 1:3
    [~, pvals(c)] = ttest(n(:,c), nCNO(:,c)); % paired t-test
end
pvals_bonf = pvals * length(pvals);%3 observation
pvals_bonf(pvals_bonf > 1) = 1 

%% Sleep Onset

%Vehicle
for s=1:numel(IntervalVeh)
Idx=find(IntervalVeh{s}(:,3)>30,1,'first');%determine first interval that last at least 30sec
if isempty(Idx)
SleepOnsetVeh(s,:)=[MaxTime*60,MaxTime*60];
else
SleepOnsetVeh(s,:)=[IntervalVeh{s}(1,1),IntervalVeh{s}(Idx,1)];
end
end
%CNO
for s=1:numel(IntervalDRUG)
Idx=find(IntervalDRUG{s}(:,3)>30,1,'first');%determine first interval that last at least 30sec
SleepOnsetCNO(s,:)=[IntervalDRUG{s}(1,1),IntervalDRUG{s}(Idx,1)];
end

COLORCODE{1}=colorcodeV;
COLORCODE{2}=colorcodeCNO;
ytitle='Sleep Onset (min)';
conditions{1}='Veh';
conditions{2}='CNO';
STATS='paired';

%% to compile data per mouse
    if Data_per_mouse ==1
    MMM=[];m=0;JJJ=[];
    Nmice=size(SleepOnsetVeh,1)/2;
        for M=1:Nmice
            MM=[];JJ=[];
            MM=nanmean(SleepOnsetVeh(1+m:2+m,1),1);
            MM(2)=nanmean(SleepOnsetVeh(1+m:2+m,2),1);
            JJ=nanmean(SleepOnsetCNO(1+m:2+m,1),1);
            JJ(2)=nanmean(SleepOnsetCNO(1+m:2+m,2),1);
            MMM=[MMM;MM];JJJ=[JJJ;JJ];
            m=m+2;
        end
        SleepOnsetVeh=MMM;
        SleepOnsetCNO=JJJ;
    end



% bar graph representation
bargraph=0;
if bargraph==1
    
figure
subplot(1,2,1)
[data1,hb1,datamean1,dataSEM1,pval1]=BarSEMgraph2([SleepOnsetVeh(:,2)./60,SleepOnsetCNO(:,2)./60], conditions, ytitle, COLORCODE,STATS);
title('min sleep 30sec','FontSize', 10)
ylim([0 100])
subplot(1,2,2)
[data2,hb2,datamean2,dataSEM2,pval2]=BarSEMgraph2([SleepOnsetVeh(:,1)./60,SleepOnsetCNO(:,1)./60], conditions, ytitle, COLORCODE,STATS);
ylim([0 100])
end

% Represent box plot with lines before after
data1=[SleepOnsetVeh(:,1)./60,SleepOnsetCNO(:,1)./60]; %transform in min - 30 sec minimum of sleep
data2=[SleepOnsetVeh(:,2)./60,SleepOnsetCNO(:,2)./60]; %transform in min


figure
subplot(1,2,1)
boxplot(data1,'labels',conditions);
hold on
color = ['r','k']; 
h = findobj(gca,'Tag','Box'); 
h2 = findobj(gca,'Tag','Median'); h3=findobj('LineStyle','--'); set(h3, 'LineStyle','-');
for j=1:length(h) 
patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.2);
h2(j).Color=color(j);
h2(j).LineWidth=2;
end 

hold on
l = length(data1(:,1));
V1 = ones(l,1)*1.25;
V2 = ones(l,1)*1.75;

for ii = 1:l
    plot([V1(ii),V2(ii)],[data1(ii,1),data1(ii,2)],'k')
end

box off
ylabel('Sleep onset (min)')
hold on 
STAT_paired_unpaired(data1,'paired');
title('min sleep 30sec','FontSize', 8)

subplot(1,2,2)
boxplot(data2,'labels',conditions);
hold on
color = ['r','k']; 
h = findobj(gca,'Tag','Box'); 
h2 = findobj(gca,'Tag','Median'); h3=findobj('LineStyle','--'); set(h3, 'LineStyle','-');
for j=1:length(h) 
patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.2);
h2(j).Color=color(j);
h2(j).LineWidth=2;
end 

hold on

for ii = 1:l
    plot([V1(ii),V2(ii)],[data2(ii,1),data2(ii,2)],'k')
end

box off
ylabel('Sleep onset (min)')
hold on 
[pval,stats]=STAT_paired_unpaired(data2,'paired');

titleFIG=['SleepOnsets_',TitleFig];
 % cd('E:\Ephy\homeCage\CNO_0.5mgkg')
 % savefig(titleFIG)

 figure
[data,datamean,dataSEM,pval]=MeanSEMgraph(data2,{'Veh','CNO'},'Sleep onset (min)',{'k','r'},'paired');
 %% # sleep bouts and ratio
% close all
 conditions={'Veh','CNO'};
 for sb=1:numel(IntervalVeh)
    NsleepBoutsVeh(sb,1)=numel(IntervalVeh{sb});
    TotTimeVeh(sb,1)=sum(IntervalVeh{sb}(:,3));
    MeanBoutsVeh(sb,1)=mean(IntervalVeh{sb}(:,3));
end
for sb=1:numel(IntervalDRUG)
    NsleepBoutsDRUG(sb,1)=numel(IntervalDRUG{sb});
    TotTimeDRUG(sb,1)=sum(IntervalDRUG{sb}(:,3));
    MeanBoutsDRUG(sb,1)=mean(IntervalDRUG{sb}(:,3));
end 

    if Data_per_mouse ==1
    MMM=[];m=0;JJJ=[];KKK=[];LLL=[];NNN=[];OOO=[];
    Nmice=size(IntervalVeh,2)/2;
        for M=1:Nmice
            MM=[];JJ=[];KK=[];LL=[];NN=[];OO=[];
            MM=nanmean(NsleepBoutsVeh(1+m:2+m,1),1);
            JJ=nanmean(NsleepBoutsDRUG(1+m:2+m,1),1);
            KK=nanmean(TotTimeVeh(1+m:2+m,1),1);
            LL=nanmean(TotTimeDRUG(1+m:2+m,1),1);
            NN=nanmean(MeanBoutsVeh(1+m:2+m,1),1);
            OO=nanmean(MeanBoutsDRUG(1+m:2+m,1),1);
            MMM=[MMM;MM];JJJ=[JJJ;JJ];KKK=[KKK;KK];LLL=[LLL;LL];NNN=[NNN;NN];OOO=[OOO;OO];
            m=m+2;
        end
        NsleepBoutsVeh=MMM;TotTimeVeh=KKK;MeanBoutsVeh=NNN;
        NsleepBoutsDRUG=JJJ;TotTimeDRUG=LLL;MeanBoutsDRUG=OOO;
    end

% Represent box plot with lines before after
data1=[NsleepBoutsVeh,NsleepBoutsDRUG]; 
data3=[TotTimeVeh./NsleepBoutsVeh,TotTimeDRUG./NsleepBoutsDRUG];%transform mean duration of bouts
data2=[MeanBoutsVeh,MeanBoutsDRUG];

figure
subplot(1,2,1)
ytitle='# bouts';
[~,datamean,dataSEM,pval]=MeanSEMgraph(data1, conditions, ytitle, {'k','r'},'paired');

subplot(1,2,2)
ytitle='Bouts duration';
[~,datamean,dataSEM,pval]=MeanSEMgraph(data2, conditions, ytitle, {'k','r'},'paired');


% %BOX PLOT Graphs:
% figure
% subplot(1,2,1)
% boxplot(data1,'labels',conditions);
% hold on
% color = ['r','k']; 
% h = findobj(gca,'Tag','Box'); 
% h2 = findobj(gca,'Tag','Median'); h3=findobj('LineStyle','--'); set(h3, 'LineStyle','-');
% for j=1:length(h) 
% patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.2);
% h2(j).Color=color(j);
% h2(j).LineWidth=2;
% end 
% 
% hold on
% l = length(data1(:,1));
% V1 = ones(l,1)*1.25;
% V2 = ones(l,1)*1.75;
% 
% for ii = 1:l
%     plot([V1(ii),V2(ii)],[data1(ii,1),data1(ii,2)],'k')
% end
% 
% box off
% ylabel(['# ',STATE,' bouts'])
% hold on 
% [pval]=STAT_paired_unpaired(data1,'paired');
% 
% subplot(1,2,2)
% boxplot(data2,'labels',conditions);
% hold on
% color = ['r','k']; 
% h = findobj(gca,'Tag','Box'); 
% h2 = findobj(gca,'Tag','Median'); h3=findobj('LineStyle','--'); set(h3, 'LineStyle','-');
% for j=1:length(h) 
% patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.2);
% h2(j).Color=color(j);
% h2(j).LineWidth=2;
% end 
% 
% hold on
% l = length(data2(:,1));
% V1 = ones(l,1)*1.25;
% V2 = ones(l,1)*1.75;
% 
% for ii = 1:l
%     plot([V1(ii),V2(ii)],[data2(ii,1),data2(ii,2)],'k')
% end
% 
% box off
% ylabel(['Mean duration of bouts in ',STATE,' (sec)'])
% hold on 
% [pval]=STAT_paired_unpaired(data2,'paired');
% 
% 
% figure
% subplot(1,2,1)
% boxplot(data3,'labels',conditions);
% hold on
% color = ['r','k']; 
% h = findobj(gca,'Tag','Box'); 
% h2 = findobj(gca,'Tag','Median'); h3=findobj('LineStyle','--'); set(h3, 'LineStyle','-');
% for j=1:length(h) 
% patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.2);
% h2(j).Color=color(j);
% h2(j).LineWidth=2;
% end 
% 
% hold on
% l = length(data3(:,1));
% V1 = ones(l,1)*1.25;
% V2 = ones(l,1)*1.75;
% 
% for ii = 1:l
%     plot([V1(ii),V2(ii)],[data3(ii,1),data3(ii,2)],'k')
% end
% 
% box off
% ylabel(['Ratio #bouts/duration in ',STATE,' (sec)'])
% hold on 
% [pval]=STAT_paired_unpaired(data3,'paired');

 %% Proportion of sleep
%see StateDistribution_Dreadd.m
end