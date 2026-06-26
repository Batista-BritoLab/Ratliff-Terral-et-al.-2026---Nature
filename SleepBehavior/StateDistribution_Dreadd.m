
function StateDistribution_Dreadd(PathlistVeh,PathlistCNO)
%Geoffrey 3/31/2026: adapted from Ratliff*,Terral* 2026
binSize=600;%in sec
% MaxTime=7200;
MaxTime=3600*4;% CHANGE MAX DURATION OF THE SESSION ACCORDINGLY
Data_per_mouse =1;

v=0;j=0;SleepScoredVeh=[];SleepScoredDRUG=[];
Pathlist=[PathlistVeh,PathlistCNO];

for s=1:numel(Pathlist)
    cd(Pathlist{s})

basename = bz_BasenameFromBasepath(pwd);
load(basename+".SleepState.states.mat")
load([basename + ".sessionInfo.mat"])

    if strcmp(sessionInfo.Treatment,'Vehicle')
        v=v+1;
       SleepScoredVeh{v}=SleepState.idx.states;
    elseif strcmp(sessionInfo.Treatment,'J60')||strcmp(sessionInfo.Treatment,'CNO')
        j=j+1;
        SleepScoredDRUG{j}=SleepState.idx.states;
    end

end

colorcodeV='k';
colorcodeCNO='r';
% NAMES=SleepState.idx.statenames;

%List of intervals of binSize until total time MaxTime
BIN=[(0:binSize:(MaxTime-binSize))',(binSize:binSize:MaxTime)'];
BIN(:,1)=BIN(:,1)+1;

%Veh
for j=1:size(SleepScoredVeh,2)
for bb=1:size(BIN,1)
%Wake proportion
WakeVeh(bb,j)=size(find(SleepScoredVeh{j}(BIN(bb,1):BIN(bb,2))==1),1)./binSize;
%Nrem proportion
NremVeh(bb,j)=size(find(SleepScoredVeh{j}(BIN(bb,1):BIN(bb,2))==3),1)./binSize;
%Rem proportion    
RemVeh(bb,j)=size(find(SleepScoredVeh{j}(BIN(bb,1):BIN(bb,2))==5),1)./binSize;
end
end

%CNO
for j=1:size(SleepScoredDRUG,2)
for bb=1:size(BIN,1)
%Wake proportion
WakeCNO(bb,j)=size(find(SleepScoredDRUG{j}(BIN(bb,1):BIN(bb,2))==1),1)./binSize;
%Nrem proportion
NremCNO(bb,j)=size(find(SleepScoredDRUG{j}(BIN(bb,1):BIN(bb,2))==3),1)./binSize;
%Rem proportion    
RemCNO(bb,j)=size(find(SleepScoredDRUG{j}(BIN(bb,1):BIN(bb,2))==5),1)./binSize;
end
end

%% to compile data per mouse
for i=1:6
    if i==1; n=WakeVeh; elseif i==2; n=NremVeh; elseif i==3; n=RemVeh; elseif i==4; n=WakeCNO; elseif i==5; n=NremCNO; elseif i==6; n=RemCNO; end
    if Data_per_mouse ==1
    MMM=[];m=0;
    Nmice=size(SleepScoredDRUG,2)/2;
        for M=1:Nmice
            MM=[];
            MM=nanmean(n(:,1+m:2+m),2);
            MMM=[MMM,MM];
            m=m+2;
        end
    if i==1; WakeVeh=MMM; elseif i==2; NremVeh=MMM; elseif i==3; RemVeh=MMM; elseif i==4; WakeCNO=MMM; elseif i==5; NremCNO=MMM; elseif i==6; RemCNO=MMM; end

    end
end

%Figure
MEANwakeVeh=mean(WakeVeh,2);
SEMwakeVeh=std(WakeVeh,[],2)./sqrt(size(WakeVeh,2)); 
MEANwakeCNO=mean(WakeCNO,2);
SEMwakeCNO=std(WakeCNO,[],2)./sqrt(size(WakeCNO,2)); 
MEANnremVeh=mean(NremVeh,2);
SEMnremVeh=std(NremVeh,[],2)./sqrt(size(NremVeh,2)); 
MEANnremCNO=mean(NremCNO,2);
SEMnremCNO=std(NremCNO,[],2)./sqrt(size(NremCNO,2)); 
MEANremVeh=mean(RemVeh,2);
SEMremVeh=std(RemVeh,[],2)./sqrt(size(RemVeh,2)); 
MEANremCNO=mean(RemCNO,2);
SEMremCNO=std(RemCNO,[],2)./sqrt(size(RemCNO,2)); 



XLABEL=(BIN(1:end,2)/60);
XLABEL=num2str(XLABEL(2:2:end));

figure
for F=1:3
subplot(1,2,1)
    if F==1%wake
        MeandataVeh=MEANwakeVeh;SEMdataVeh=SEMwakeVeh;MeandataCNO=MEANwakeCNO;SEMdataCNO=SEMwakeCNO;COLOR='k';
    elseif F==2%nrem
        MeandataVeh=MEANnremVeh;SEMdataVeh=SEMnremVeh;MeandataCNO=MEANnremCNO;SEMdataCNO=SEMnremCNO;COLOR='b';
    elseif F==3%rem
        MeandataVeh=MEANremVeh;SEMdataVeh=SEMremVeh;MeandataCNO=MEANremCNO;SEMdataCNO=SEMremCNO;COLOR='r';
    end
    
hb1=shadedErrorBar([1:size(BIN,1)],smooth(MeandataVeh,0.1,'loess'),SEMdataVeh);
hb1.mainLine.Color=COLOR;
hb1.patch.FaceColor=COLOR;
hb1.patch.FaceAlpha=0.5;
set(hb1.edge,'LineWidth',2,'LineStyle','none');

hold on
subplot(1,2,2)
hb2=shadedErrorBar([1:size(BIN,1)],smooth(MeandataCNO,0.1,'loess'),SEMdataCNO);
hb2.mainLine.Color=COLOR;
hb2.patch.FaceColor=COLOR;
hb2.patch.FaceAlpha=0.5;
set(hb2.edge,'LineWidth',2,'LineStyle','none');
end
subplot(1,2,1)
title('Vehicle')
ylabel('% time in state')
xticks(2:2:size(XLABEL,1)*2)
xticklabels(XLABEL);

subplot(1,2,2)
title('CNO')
legend('WAKE','NREM','REM')
xticks(2:2:size(XLABEL,1)*2)
xticklabels(XLABEL);

end