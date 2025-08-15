function [outcome, curLickEvents]=Online_LickEvents(currentTrial,StateToZero)
%[outcome, curLickEvents]=CurrentTrialEvents(BpodSystem, trialsMatrix, currentTrial, currentTrialType, time)
%
%This function extracts the outcome (absence or presence of neverlickedstate) and the licking events
%to update the trials and licks plots, respectively (see associated functions). 
%The timestamp of lickEvents output is normalized to the timing of the event 
%specified by the input argument "type" (cue or reward),
%
%Output arguments can be used as an input argument for Online_LickPlot function.
%
%function written by Quentin for CuedReinforcers bpod protocol

global BpodSystem
% %% Extract the outcome to update the trialsplot function
% if trialsMatrix(currentTrialType,6)==1;
%     if isnan(BpodSystem.Data.RawEvents.Trial{1,currentTrial}.States.NeverLickedState(1,1))==1
%         outcome='g'; %hit : did not go to neverlickedstate
%     else
%         outcome='r'; %missed :went to neverlickedstate
%     end
% else    outcome='g'; %uncued reward
% end

%% Extract the lick events from the BpodSystem structure
curLickEvents=NaN;  %if no lick, random number
try
    LickEventsRaw=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.Events.Port1In;
    TimeForZero=BpodSystem.Data.RawEvents.Trial{1, currentTrial}.States.(StateToZero)(1,1);      
    curLickEvents=LickEventsRaw-TimeForZero;
end
%% Outcome : green if the animal has collected the reward / showed anticipatory licks
if sum(curLickEvents>-1 & curLickEvents<2)>1
    outcome='g';
else
    outcome='r';
end

end

