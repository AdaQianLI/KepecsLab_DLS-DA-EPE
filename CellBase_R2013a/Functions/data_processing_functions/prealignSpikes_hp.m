function  prealignSpikes_hp(mycellids,varargin)
%
%  PREALIGNSPIKES       Aligns spikes to events and calculates rates for
%                       each epoch.
%
%
%  prealignSpikes(cellid(s),{'default'|'defineEventsEpochs'})
%
%  If the content of cellids is 'all' then the function will run on all the
%  cells in CellBase.
%
%  You should create your own event/epoch def. eg.: @defineEventsEpochs_adam1106
%
% STIMES data format:
% ST.CELLID
% ST.EVENTS
% ST.EVENT_STIMES
% ST.EVENT_WINDOWS
% ST.EPOCHS
% ST.EPOCH_RATES

% AK 11/06

default_args={...
    'ifsave',      0;...
    'ifappend',    1;...
    'FUNdefineEventsEpochs', @defineEventsEpochs_default;...
    'filetype',    'behav';... % 'event' 'stim' This determines what type of events we are prealigning to
    };
[g, error] = parse_args(default_args,varargin{:});

% if (nargin == 1) || strcmpi(varargin{1},'default')
%     FUNdefineEventsEpochs = @defineEventsEpochs_default;
%     disp('PrealignSpikes is using default event and epoch definitions.')
% else
%     FUNdefineEventsEpochs = varargin{1};
% end

% This is for addnewcells, so that we can use the latest version of events
% and epochs to prealignspikes for the newcells.
if ~isfield(g,'events') %&& ~isfield(g,'epochs'),
   [events,epochs] = feval(g.FUNdefineEventsEpochs);
else
    events=g.events;
    %epochs=g.epochs;
end
switch g.filetype
    case {'behav','event'}
        FNAMESTR='EVENTSPIKES';
        EVFILE='TrialEvent';
    case 'stim'
        FNAMESTR='STIMSPIKES';
        EVFILE='StimEvent';
end
% FNAMESTR = 'STIMES';
%FNAMESTR = 'EVENTSPIKES';

if strcmpi(mycellids,'all')
      mycellids = listtag('cells');
end
    
for i=1:length(mycellids)
   cellid = mycellids(i);
   EV=loadcb(cellid,EVFILE);
   % now we need to pass the Event file to extractEventSpikes
   [event_stimes,event_windows] = extractEventSpikes_hp(cellid,events,EV);
   %[epoch_rates,dur,count] = extractEpochRates(event_stimes,event_windows,events,epochs);

   fname = cellid2fnames(cellid,FNAMESTR);
   if g.ifappend==0
       try
           %save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates','-append');
           save(fname,'cellid','events','event_stimes','event_windows'); % no need for epochs
       catch
           disp('No file to append')
           %save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates');
           save(fname,'cellid','events','event_stimes','event_windows'); % no need for epochs
       end
   else
       %save(fname,'cellid','events','event_stimes','event_windows','epochs','epoch_rates');
       save(fname,'cellid','events','event_stimes','event_windows','-append'); % no need for epochs
   end
end