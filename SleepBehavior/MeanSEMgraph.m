
function [data,datamean,dataSEM,pval]=MeanSEMgraph(data, conditions, ytitle, COLORCODE,STATS)
% STATS='paired';%tell if paied or unpaired - if empty no running stat
% STATS='unpaired';%1 if you want to display stat - by default
% Find color by default
% % oldFolder=pwd;
% % cd('D:\')
% % load('ColorCode.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lengthLine=0.1;

if ~iscell(data)
[row,col]=find(data==inf);
if ~isempty(row)
for i=1:size(row,1)
data(row(i),col(i))=NaN; %in case Inf
end 
end
data_mean=nanmean(data,1);

NanVal=size(data,1)-sum(isnan(data),1);%how many units not to considered
data_SEM = std(data,[],1,'omitnan')./sqrt(NanVal); 
else %data need to be in columns
    data_mean=[];data_SEM=[];
    for C=1:size(data,2)
[row,col]=find(data{:,C}==inf);
if ~isempty(row)
for i=1:size(row,1)
    data{:,C}(row(i),col(i))=NaN; %in case Inf
end 
end
% data_mean=[data_mean,nanmedian(data{:,C},'all')];
data_mean=[data_mean,nanmean(data{:,C},'all')];

NanVal=size(data{:,C},1)-sum(isnan(data{:,C}),1);%how many units not to considered
data_SEM = [data_SEM,std(data{:,C},[],1,'omitnan')/sqrt(NanVal)]; 
    end

end

X=1:size(data_mean,2);


% figure,
% axisLimits = axis;
% lengthInDataUnits = diff(axisLimits(1:2)) * lengthLine;

if strcmp(STATS,'paired') && size(data_mean,2)==2
p=plot([X(1),X(2)],[data_mean(1),data_mean(2)],'LineWidth',4);
p.Color=[0.5 0.5 0.5];
end

hold on;

for i= 1:size(data_mean,2)

% for i= 1:length(data_mean) %bug 24/09/21 with cell array in data
% hb(i)=bar (i,data_mean(i));%here to change to get bar graph
% hb(i)=plot(i,data_mean(i),'Marker','o','Color',COLORCODE{i});
% hb(i) = line([X(i)-lengthInDataUnits/2, X(i)+lengthInDataUnits/2], [data_mean(i), data_mean(i)], 'Color', COLORCODE{i}, 'LineWidth', 3);
% 
% hold on
% % set(hb(i),'FaceColor',COLORCODE{i},'FaceAlpha',0.2);%here to change to get bar graph
% set(hb(i),'MarkerFaceColor',COLORCODE{i},'MarkerSize', 8,'MarkerEdgeColor', 'none');

% Plot error bars for the left point with black color
er = errorbar(X(i), data_mean(i), data_SEM(i), 'LineStyle', 'none', 'CapSize', 16, 'LineWidth', 4, 'Color', COLORCODE{i});
% %To show mean with a bar line
% er2 = errorbar(X(i), data_mean(i), 0,'LineStyle', 'none', 'CapSize', 14, 'LineWidth', 3, 'Color', COLORCODE{i});
% set (gca, 'XTickLabel',conditions{i}, 'Fontsize', 14)
% set (gca, 'XTickLabel',conditions{i})

end

hold on

% er= errorbar(X,data_mean,data_SEM);
% er.LineStyle= 'none';
% er.CapSize=8;
% er.LineWidth=2;


% set (gca,'Xtick', X, 'XTickLabel',conditions, 'Fontsize', 14)
 set (gca,'Xtick', X, 'XTickLabel',conditions,'Fontsize', 18)
xtickangle(45) %%rotate tick of 45d
ylabel(ytitle,'Fontsize',18)
box off

datamean=data_mean;
dataSEM=data_SEM;
% % %return to previous folder
% % cd(oldFolder)

% % if scatter plot as well on top of the graph
hold on

    if strcmp(STATS,'paired')
    hold on
    if iscell(data)
        data=[data{1} data{2}];%if paired it cannot be cell
    end
    l = length(data(:,1));
    V1 = ones(l,1)*(1+lengthLine+lengthLine/10);
    V2 = ones(l,1)*(2-lengthLine-lengthLine/10);

    for j = 1:l
        % if ~iscell(data)
    pp=plot([V1(j),V2(j)],[data(j,1),data(j,2)],'LineWidth',0.15);
    pp.Color=[0.5 0.5 0.5];
    %     else
    % plot([V1(j),V2(j)],[data{j,1},data{j,2}],'k','LineWidth',0.5);
    %     end
    end


    hold on 
    [pval,~]=STAT_paired_unpaired([data],'paired');
    else

    for i= 1:size(data_mean,2)
    if ~iscell(data)
    scatter(ones(size(data(:,i),1),1)*i,data(:,i),10,'MarkerFaceColor',COLORCODE{i},'MarkerEdgeColor','none')
    else
    scatter(ones(size(data{i},1),1)*i,data{i},10,'MarkerFaceColor',COLORCODE{i},'MarkerEdgeColor','none')
    end
    end
    hold on 
    [pval,~]=STAT_paired_unpaired([data],'unpaired');  

    end

    box off

% %% Stats
% %Stats
% if nargin<5
%     STATS=[];
% end
% 
% if ~isempty(STATS)
% if size(data_mean,2)>2 %in case 3 column - specify in conditions instead the text of injection
% 
% A=[];
% for u=1:3  
%     if iscell(data)
%     v=ones(size(data{u}))*u; %work if in cells
%     A=[A;[data{u},v]];
%     else
%     v=ones(size(data(:,u)))*u; 
%     A=[A;[data(:,u),v]];
%     end
% 
% end
% 
% 
% [pval,~,stats] = kruskalwallis(A(:,1),A(:,2),'off');
% if pval<=0.05
%    c = multcompare(stats ,'displayopt','off','ctype','bonferroni');
% 
% if c(1,6)>0.05; c(1,6)=nan;end
% if c(2,6)>0.05; c(2,6)=nan;end
% if c(3,6)>0.05; c(3,6)=nan;end
% H=sigstar({[c(1,(1:2))],[c(2,(1:2))],[c(3,(1:2))]},[c(1,6),c(2,6),c(3,6)]);
% set(H(:,2),'FontSize', 15);
% elseif pval>0.05 ;p=nan; 
%     H=sigstar({[1,3]},p);set(H(:,2),'FontSize', 15);
% end
% 
% else
%     if strcmp(STATS,'unpaired')
%         if iscell(data)
%         [pval,~,stats]=ranksum(data{1},data{2});      
%         else
%         [pval,~,stats]=ranksum(data(:,1),data(:,2));
%         end
%     elseif strcmp(STATS,'paired')    
%         if iscell(data)
%         [pval,~,stats]=signrank(data{1},data{2});      
%         else
%         [pval,~,stats]=signrank(data(:,1),data(:,2));
%         end
%     end
%         if pval<=0.05
%         H=sigstar({[1 2]},pval);
%         else
%             p=nan;
%         H=sigstar({[1 2]},p);
%         end
%         set(H(:,2),'FontSize', 15);
% end
% end
xlim([0 size(data_mean,2)+1])
end