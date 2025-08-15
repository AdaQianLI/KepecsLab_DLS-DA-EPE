function [outcome, curLickEvents]=Online_LickEvents_8class(currentTrial,StateToZero)
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

TimeForZero=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.States.(StateToZero)(1,1);    
    try
       LickEventsRaw=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.Events.Port1In; 
    catch
       LickEventsRaw=NaN;  %if no lick, random number
    end
       curLickEvents=LickEventsRaw-TimeForZero;
%% Outcome : green if the animal has collected the reward / showed anticipatory licks
  try
       MotionEventsRaw=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.States.WaterDeliver(1,1);
  catch
       MotionEventsRaw=NaN;
  end
      if ~isnan(MotionEventsRaw)
          outcome='g';
      else 
          outcome='r';
      end
end

