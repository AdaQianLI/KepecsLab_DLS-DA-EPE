function TrialType=Online_TrialType_8class(currentTrial,seq)
% ONLINE_TRIALTYPE_8CLASS
% Convert the programmed trial identity plus the actual outcome into one of
% eight result classes used by the online figures.
global BpodSystem


%% Reward delivery determines whether the current trial is labeled success or fail

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

