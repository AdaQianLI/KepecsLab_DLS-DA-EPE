function [outcome, curLickEvents]=Online_LickEvents_8class(currentTrial,StateToZero)
% ONLINE_LICKEVENTS_8CLASS
% Extract lick timestamps and trial outcome from the current trial so the
% online lick raster and trial-outcome plot can be updated.
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
%% Lick timestamps aligned to the chosen event

TimeForZero=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.States.(StateToZero)(1,1);    
    try
       LickEventsRaw=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.Events.Port1In; 
    catch
       LickEventsRaw=NaN;  % No licks recorded on this trial.
    end
       curLickEvents=LickEventsRaw-TimeForZero;
%% Outcome flag for the online trial plot
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

