function [ranova_results,BonfStat]=figureCalciumActivity_inStates(CALCIUM, ylabeltitle)
% figure
% plot(Calcium(s).time,deltaF,'k')

%With REM
figure
[data2,hb,datamean,dataSEM,pval]=BarSEMgraph2(CALCIUM(:,1:4),{'SWS','Quiet wake','Movement','REM'}, ylabeltitle, {'k','k','k','k'},[]);
hold on
hb=scatter([1,2,3,4],CALCIUM,'k','filled');
hold on
for ii = 1:size(CALCIUM,1)
plot([1,2,3,4],CALCIUM(ii,:),'k')
end
xlim([0 5])
box off



%% ANOVA:

% Define the repeated-measures model
if anynan(CALCIUM(:,4)) && ~all(isnan(CALCIUM(:,4)))
    [ranova_results,BonfStat]=fitlmeCalcium(CALCIUM); %mixed-effects model - if only few missing data NaN
    return;
end

if all(isnan(CALCIUM(:,4)), 1) 
    ndata=[1:3]';
    ANOVAterm='SWS-Move ~ 1';
    tbl = array2table(CALCIUM(:,1:3), 'VariableNames', {'SWS', 'QW', 'Move'});

elseif all(~isnan(CALCIUM(:,4)), 1) %if all nan because no REM or a
    ndata=[1:4]';
    ANOVAterm='SWS-REM ~ 1';
    tbl = array2table(CALCIUM(:,1:4), 'VariableNames', {'SWS', 'QW', 'Move','REM'});
end

pairs = GetPairsFromList_LR(ndata);

rm = fitrm(tbl, ANOVAterm, 'WithinDesign', ndata);
% Run repeated-measures ANOVA
ranova_results = ranova(rm);
% Display results
% disp(ranova_results);
BonfStat = [];
if ranova_results.pValue(1) <= 0.05

% Column pairs to compare (1 vs 2, 2 vs 3, 1 vs 3)
% pairs = [1 2; 2 3; 1 3; 1 4; 2 4; 3 4];
for i = 1:size(pairs,1)
    condA = pairs(i,1);
    condB = pairs(i,2);
    %paired t-test
    [~, p] = ttest(CALCIUM(:,condA), CALCIUM(:,condB));
    BonfStat = [BonfStat; condA, condB, p];
end

% Bonferroni Correction (divide p-values by number of comparisons)
numComparisons = size(pairs,1);
BonfStat(:,3) = min(BonfStat(:,3) * numComparisons, 1); %here give multiple comparision for specific pairs 
end

end


function [ranova_results,BonfStat]=fitlmeCalcium(CALCIUM) %mixed-effects model - if only few missing data NaN
pairs = [1 2; 2 3; 1 3; 1 4; 2 4; 3 4];
BonfStat = [];
[nSubjects, nCond] = size(CALCIUM);
Response = reshape(CALCIUM(:,1:nCond)', [], 1);   
Subject  = repelem((1:nSubjects)', nCond);       
conds = {'SWS','QW','Move','REM'};
Condition = repmat(conds, 1, nSubjects)';      

tbl_long = table(Response, categorical(Subject), categorical(Condition), ...
                 'VariableNames', {'Response','Subject','Condition'});
tbl_long.Condition = reordercats(tbl_long.Condition, conds);%here change to get SWS as "baseline" Intercept
lme = fitlme(tbl_long, 'Response ~ Condition + (1|Subject)');

ranova_results = anova(lme, 'DFMethod','Satterthwaite');
if ranova_results.pValue(2) <= 0.05
for i = 1:size(pairs,1)
    condA = pairs(i,1);
    condB = pairs(i,2);

    % Build contrast: coefTest does the paired comparison
    names = lme.CoefficientNames;
    L = zeros(1,numel(names));

    if condA == 1 % SWS baseline
        L(strcmp(names, ['Condition_' conds{condB}])) = 1;
    elseif condB == 1
        L(strcmp(names, ['Condition_' conds{condA}])) = 1;
        L = -L;
    else
        L(strcmp(names, ['Condition_' conds{condA}])) = 1;
        L(strcmp(names, ['Condition_' conds{condB}])) = -1;
    end

    [p,~,~,~] = coefTest(lme,L);
    BonfStat = [BonfStat; condA, condB, p];
end
% Bonferroni correction
BonfStat(:,3) = min(BonfStat(:,3) * size(pairs,1), 1);
else
    BonfStat=[];

end
end