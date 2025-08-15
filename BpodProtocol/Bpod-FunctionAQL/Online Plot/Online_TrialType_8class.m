function TrialType=Online_TrialType_8class(currentTrial,seq)
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

try
CurTrialType=BpodSystem.Data.RawEvents.Trial{1,currentTrial}.States.WaterDeliver(1,1);
catch
    CurTrialType=[];
end
switch seq
    case 1
           if ~isnan(CurTrialType)
            TrialType=1;
           else
            TrialType=2;
           end
           
    case 2
           if ~isnan(CurTrialType)
            TrialType=3;
           else
            TrialType=4;
           end
    case 3
           if ~isnan(CurTrialType)
            TrialType=5;
           else 
            TrialType=6;
        end
    case 4
           if ~isnan(CurTrialType)
            TrialType=7;
           else
            TrialType=8;
          end
end
end

