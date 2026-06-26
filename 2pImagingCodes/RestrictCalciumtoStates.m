
function CALCIUM=RestrictCalciumtoStates(Calcium,statesIntervals,IntanOnset)

%IntanOnset= treadmill.time(1);

RestrictedTPS=[];
tot_frame=Calcium.frameTransition;
apprFramerate=mean(diff(Calcium.time));
BIN=round(600/apprFramerate);%1min window

for ii=1:size(Calcium.signal,1)
    XmaxquantMIN = movquant(Calcium.rawF(ii,1:tot_frame), 0.1, BIN, [2], 'omitnan', 'truncate');
    XmaxquantMAX = movquant(Calcium.rawF(ii,1:tot_frame), 0.9, BIN, [2], 'omitnan', 'truncate');

    deltaF=(Calcium.rawF(ii,1:tot_frame) - XmaxquantMIN)./(XmaxquantMIN); %Same normalization as Jacob - Manuscript 2024
    
RestrictedTPS=InIntervals(Calcium.time(1:tot_frame),statesIntervals.SWS+IntanOnset);
CALCIUM(ii,1)=nanmean(deltaF(RestrictedTPS));
RestrictedTPS=InIntervals(Calcium.time(1:tot_frame),statesIntervals.QuietWAKE+IntanOnset);
CALCIUM(ii,2)=nanmean(deltaF(RestrictedTPS));
RestrictedTPS=InIntervals(Calcium.time(1:tot_frame),statesIntervals.Movement+IntanOnset);
CALCIUM(ii,3)=nanmean(deltaF(RestrictedTPS));

 if ~isempty(statesIntervals.REM)
    RestrictedTPS=InIntervals(Calcium.time(1:tot_frame),statesIntervals.REM+IntanOnset);
    CALCIUM(ii,4)=nanmean(deltaF(RestrictedTPS));    
 else
    CALCIUM(ii,4)=NaN;    
 end

% CALCIUM(ii,5)=Calcium(s).Depth(ii,1);%depth or instead type1/type2 to
% add?
end

end