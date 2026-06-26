function Sync_Imaging_down_Ratliff_Terral2026

%See also original script extractPeriEventSignal
fig=1;

%For shuffle:
nShuffle = 100;

%Default values
% if nargin<1
    basepath=pwd;
% end

% basepath=pwd;
basename = bz_BasenameFromBasepath(basepath);
load(basename + ".outputStruct.mat");%
load(basename + ".SlowWaves.events.mat");%
load(basename + ".SleepState.states.mat");%
NREM=SleepState.ints.NREMstate;
% get exact value for intan (outputStruct is downsample 100Hz)
load(basename + ".ImgStruct.mat");%
timeIntan=imgStruct.AlignedData.treadmillData.time';
% fs=round(1/(timeIntan(2)-timeIntan(1)));

SOtps=SlowWaves.ints.DOWN;
SUtps=SlowWaves.ints.UP;
SUtps=SUtps+timeIntan(1);
SUtps(SUtps(:,1)<0,:)=[];

% smoothWin = 20; % frames
smoothWin=1;
% smoothWin = fs*1.5; % 1.5 sec smoothing

% signalImag=img.rawF;
signalImag=img.Deconv(:,1:img.frameTransition);%deconvolution
% signalImag=img.signal(:,1:img.frameTransition);%DF/F
% 
SOtps=SOtps+timeIntan(1);
SOtps(SOtps(:,1)<0,:)=[];
NREM=NREM+timeIntan(1);
NREM(NREM(:,1)<0,:)=[];
NREM=IntersectIntervals(NREM,[1 img.time(end)]);

timeIntan = timeIntan(timeIntan > 0 & timeIntan <= img.time(end));
timeIntan=downsample(timeIntan,20);
fs=round(1/mean(diff(timeIntan)));
NREMidx = zeros(size(timeIntan));

%2 second exclusion; then remove 0.5sec before/after each interval to avoid problem in getting random values during NREM before/after down
nexclNREM=NREM;
nexclNREM(NREM(:,2)-NREM(:,1)<2,:)=[];%3 second exclusion: will remove 1sec 
nexclNREM=[nexclNREM(:,1)+0.5,nexclNREM(:,2)-0.5];

NREMidx=[];
for k = 1:size(nexclNREM,1)
    idx = find(timeIntan >= nexclNREM(k,1) & timeIntan <= nexclNREM(k,2));
    NREMidx = [NREMidx;idx];
end
% NREMidx(randi(length(NREMidx)));

%extrapolate imaging to intan

ImgDataStreched(1,:) = interp1(img.time, signalImag, timeIntan);

% timeIntan=timeIntan(NREMidx);%8/21/25
% % %Zscore
% ImgDataStreched=ImgDataStreched(NREMidx);

ImgDataStreched=zscore(ImgDataStreched);
%Norm to max
% ImgDataStreched=ImgDataStreched./max(ImgDataStreched);

signalImag=smooth(ImgDataStreched,smoothWin)';

% numPointsBefore = 10000; % Number of points before the interval
% numPointsOut = fs/2;%in sampling(500msec)
% numPointsInside = fs/2; % Number of points inside the interval
numPointsOut = fs;%in sampling(1sec)
numPointsInside = fs/2; % Number of points inside the interval

periEventMatrix3=[];IDX_Down_used=[];Downint=[];
SzIm=size(signalImag,2);
k=1;

tic

% [SOtps,~]=RestrictInts(SOtps,nexclNREM);%RESTRICT DOWN IN NREM

%Restrict to events long enough according to 1sec
RestSotps=zeros(size(SOtps,1),1);
for i=2:size(SOtps,1)
    a=(SOtps(i,1)-(SOtps(i-1,2)));
    if a<1
    RestSotps(i,1)=false;
    else
        RestSotps(i,1)=true;
    end
end
RestSotps(1,1)=true;
        
SOtps=SOtps(logical(RestSotps),:);

periEventMatrix_shuffle = nan(1, 2*numPointsOut+numPointsInside, nShuffle);

for i=1:size(SOtps,1)
    startT = SOtps(i, 1);
    endT = SOtps(i, 2);
    Idximg=find(timeIntan<=startT,1,'last');
    IdximgOffset=find(timeIntan<=endT,1,'last');
    % Idximg2=find(img.time<=startT,1,'last');
    % IdximgOffset2=find(img.time<=endT,1,'last');

    if ~isempty(Idximg) && ~isempty(IdximgOffset)

        if (IdximgOffset+numPointsOut)<=SzIm && (Idximg-numPointsOut)>=0
           % idxInside = find(timeIntan >= startT & timeIntan <= endT); 
           idxInside = [Idximg:IdximgOffset];
           % idxInside2 = [Idximg2:IdximgOffset2];
           resampledTime = linspace(startT, endT, numPointsInside);
           resampledData = interp1(timeIntan(idxInside), signalImag(idxInside), resampledTime, 'linear', 'extrap');

           periEventMatrix3(k,:)=[signalImag(1,Idximg-numPointsOut:Idximg-1),resampledData,signalImag(1,(IdximgOffset+1):(IdximgOffset+1)+numPointsOut-1)];        


           Downint=[Downint;startT,endT];
           IDX_Down_used=[IDX_Down_used;i];


                    for s = 1:nShuffle
                    % shift = NREMidx(randi(length(NREMidx)));
                    shift = randi([-fs fs]);                  
                    % shift = allowedShifts(randi(length(allowedShifts)));
                    % shift = SUtpsidx(randi(length(SUtpsidx)));
                    signalImag_shuff = circshift(signalImag, shift);%here should restrict img to NREM
                    % signalImag_shuff = circshift(signalImagSU, shift);
                    resampledData = interp1(timeIntan(idxInside), signalImag_shuff(idxInside), resampledTime, 'linear', 'extrap');
                    periEventMatrix_shuffle(k,:,s) = [signalImag_shuff(1,Idximg-numPointsOut:Idximg-1), resampledData, signalImag_shuff(1,(IdximgOffset+1):(IdximgOffset+1)+numPointsOut-1)];
                    end

                k=k+1;

        end

    end      

end

disp(['to run NShuffle = ', num2str(nShuffle)])
toc

mean_shuff = mean(periEventMatrix_shuffle,3);  
mean_shuff = squeeze(mean_shuff); 

% compute mean and SEM across shuffles
shuff_mean = mean(mean_shuff,1);
shuff_sem  = std(mean_shuff,0,1)/sqrt(k-1);


%% figure
if fig==1
h=figure;
% hb=shadedErrorBar(timeIntan(1:numPointsBefore+numPointsInside),nanmean(allIntervals'),std(allIntervals')./sqrt(size(allIntervals',1)));- comment 8/20/25

hb=shadedErrorBar(timeIntan(1:(2*numPointsOut+numPointsInside))-timeIntan(1),nanmean(periEventMatrix3),std(periEventMatrix3)./sqrt(size(periEventMatrix3,1)));

hb.mainLine.Color='b';
hb.patch.FaceColor='b';
hb.patch.FaceAlpha=0.5;
set(hb.edge,'LineWidth',2,'LineStyle','none');
ylabel('Zscore DF/F')
xlabel('Time (sec)')


hold on
% hb2=shadedErrorBar(timeIntan(1:(2*numPointsOut+numPointsInside)),nanmean(periEventMatrix_shuffle),std(periEventMatrix_shuffle)./sqrt(size(periEventMatrix_shuffle,1)));
hb2 = shadedErrorBar(timeIntan(1:(2*numPointsOut+numPointsInside))-timeIntan(1), shuff_mean, shuff_sem);

hb2.mainLine.Color='k';
hb2.patch.FaceColor='k';
hb2.patch.FaceAlpha=0.5;
set(hb2.edge,'LineWidth',2,'LineStyle','none');


xline(timeIntan(numPointsOut)-timeIntan(1), 'LineStyle', '--', 'LineWidth', 0.5); 
xline(timeIntan(numPointsOut+numPointsInside)-timeIntan(1), 'LineStyle', '--', 'LineWidth', 0.5); 

%Subtract data by shuffle
periEventMatrix4=periEventMatrix3-mean_shuff;
h2=figure;
% hb=shadedErrorBar(timeIntan(1:numPointsBefore+numPointsInside),nanmean(allIntervals'),std(allIntervals')./sqrt(size(allIntervals',1)));- comment 8/20/25

hb=shadedErrorBar(timeIntan(1:(2*numPointsOut+numPointsInside))-timeIntan(1),nanmean(periEventMatrix4),std(periEventMatrix4)./sqrt(size(periEventMatrix4,1)));

hb.mainLine.Color='b';
hb.patch.FaceColor='b';
hb.patch.FaceAlpha=0.5;
set(hb.edge,'LineWidth',2,'LineStyle','none');
ylabel('Change over shuffle - ZDeconvolved Ca+')
xlabel('Time (sec)')
xline(timeIntan(numPointsOut)-timeIntan(1), 'LineStyle', '--', 'LineWidth', 0.5); 
xline(timeIntan(numPointsOut+numPointsInside)-timeIntan(1), 'LineStyle', '--', 'LineWidth', 0.5); 
end
% 

graphTimeCourse.x=timeIntan(1:(2*numPointsOut+numPointsInside))-timeIntan(1);
graphTimeCourse.values=periEventMatrix3;
graphTimeCourse.chance=periEventMatrix_shuffle;


%% Calculate Average values before during after DOWN:
ImgUPonset=nanmean(periEventMatrix3(:,1:numPointsOut),2);
ImgUPoffset=nanmean(periEventMatrix3(:,(numPointsOut+numPointsInside+1):(numPointsInside+numPointsOut*2)),2);
ImgDOWN=nanmean(periEventMatrix3(:,numPointsOut+1:(numPointsOut+numPointsInside)),2);

ImgUPonset_chance=squeeze(nanmean(periEventMatrix_shuffle(:,1:numPointsOut,:),2)); 
ImgUPoffset_chance=squeeze(nanmean(periEventMatrix_shuffle(:,(numPointsOut+numPointsInside+1):(numPointsInside+numPointsOut*2),:),2));
ImgDOWN_chance=squeeze(nanmean(periEventMatrix_shuffle(:,numPointsOut+1:(numPointsOut+numPointsInside),:),2));


if fig==12
figure
pval=Boxplot_GT([ImgUPonset,ImgDOWN],{'up-on','down'},{'r','k'},'paired');
figure
pval=Boxplot_GT([ImgUPoffset,ImgDOWN],{'up-off','down'},{'r','k'},'paired');
figure
pval=Boxplot_GT([(ImgUPonset+ImgUPoffset)./2,ImgDOWN],{'up','down'},{'r','k'},'paired');

end

DurationDOWN=Downint(:,2)-Downint(:,1);
AmpDOWN=SlowWaves.SWpeakmag(IDX_Down_used,:);

%IDX_Down_used= from the raw down intervals, which one did we used for this
%analysis( keep only if signal present before and after - not cut in middle)
save (basename + ".DownSimaging1sec.mat",'ImgUPoffset','ImgUPoffset_chance','ImgUPonset','ImgUPonset_chance','ImgDOWN','ImgDOWN_chance','DurationDOWN','IDX_Down_used','graphTimeCourse','Downint','AmpDOWN','nShuffle')
% save (basename + ".DownSimaging.mat",'ImgUP','ImgUPoffset','ImgUPonset','ImgDOWN','DurationDOWN','DurationUP','graphTimeCourse','Downint','UPint','SpaceDOWN')
if fig==1
%save figure
savefig(h, fullfile(basepath, 'Imaging_DownStatesPSTH1sec.fig'));
% Save as .pdf
exportgraphics(h, fullfile(basepath, 'Imaging_DownStatesPSTH1sec.pdf'), 'ContentType', 'vector'); % vector graphics for quality
end

end
