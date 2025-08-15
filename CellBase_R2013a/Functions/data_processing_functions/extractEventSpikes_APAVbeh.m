function trial_stimes= extractEventSpikes_APAVbeh(cellid,events,win)
%
%   EXTRACTEVENTSPIKES      Extracts and aligns spikes in each trial
%                           relative to trials events
%
%   [event_stimes,event_windows] = extractEventSpikes(cellid,events)
%
%   cellid  - id of a unit in CellBase
%   events  - Nx4 cell array of events (see defineEventsEpochs)
%
%   event_stimes - 1xN cell array of 1xTrialNum cell arrays of spike times
%                  relative to the first event trigger in 'events'.
%   event_window - 1xN cell array of 1xTrialNum arrays noting the time window 
%                  of extraction for each trial.


%   AK 11/06, ZFM 10/04


disp(['Extracting event-triggered spike times for ' char(cellid)])

NUMevents     = length(events);
% event_stimes  = cell(NUMevents,1);
% event_windows = cell(NUMevents,1);

%load TrialEvents
% if nargin<3,
%     TE = loadcb(cellid,'Events');   
% end

%load spike times
stimes = loadcb(cellid,'Spikes');
%not sure why this was here -AK; n
%stimes = stimes * 1e-4; % conversion factor into seconds

% for iE=1:NUMevents
    
%     % event references 
%     TriggerEvent1 = events{iE,2}; 
%     TriggerEvent2 = events{iE,3};
%     
%     % window 
%     offset1 = events{iE,4}(1);
%     offset2 = events{iE,4}(2);
    
%     if isfield(TE,TriggerEvent1) & isfield(TE,TriggerEvent2)
        
%         NUMtrials = length(getfield(TE,TriggerEvent1));
%         
%         
%         EventTimes1 = getfield(TE,TriggerEvent1)  + offset1;
%         EventTimes2 = getfield(TE,TriggerEvent2)  + offset2;
        %EventTimes1 = getfield(TE,TriggerEvent1) + TE.TrialStart + offset1;
        %EventTimes2 = getfield(TE,TriggerEvent2) + TE.TrialStart + offset2;        
        
        % note event windows are relative to TriggerEvent1
%         event_windows{iE}(1,1:NUMtrials) = offset1;  
%         event_windows{iE}(2,:)           = (EventTimes2-EventTimes1)+offset1;
        
        
        trial_stimes = cell(1,NUMevents);
        for iT=1:NUMevents      
       if ~isnan(win{1,iT})        
            %if ~isnan(EventTimes1(iT)) & ~isnan(EventTimes2(iT))    
                temp= stimes((stimes >= win{iT}(1) & stimes < win{iT}(2)));
                trial_stimes{iT}=temp-events(iT);
       end
%         event_stimes{iE,:} = trial_stimes;
        end %iT  
%     else % some events not found
%         
%         disp(sprintf('Warning: %s missing events: %s and/or %s',char(cellid),TriggerEvent1,TriggerEvent2));
%         event_stimes{iE} = [];
%         event_windows{iE} = [NaN NaN];
% %     end
% end %iE

disp('hello')

end