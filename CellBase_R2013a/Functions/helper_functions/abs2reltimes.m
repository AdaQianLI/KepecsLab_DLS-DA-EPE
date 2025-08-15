function reltimes = abs2reltimes(cellid,varargin)
%ABS2RELTIMES   Calculate relative time windows relative to events.
%   RELTIMES = ABS2RELTIMES(CELLID,ABSTIMES,EVENT_TYPE,EVENT) calculates
%   time stamps (RELTIMES) corresponding to ABSTIMES (absolute time stamps,
%   in seconds) relative to an event (EVENT) for a given cell (CELLID).
%   EVENT_TYPE should determine whether EVENT is a field of stimulus events
%   ('stim') or trial events ('trial'). For each absolute time stamp, only
%   one relative time stamp is returned, relative to the closest event of
%   the specified type before it.
%
%   See also REL2ABSTIMES and FINDSEGS2.

%   Edit log: BH 4/25/12, 5/4/12

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'abstimes',@(s)isnumeric(s)) % time stamps to convert
addParameter(prs,'event_type','stim',@ischar)   % event type - stim or trialparse(prs,cellid,varargin{:})
addParameter(prs,'event','PulseOn_chr2',@ischar)   % reference event
parse(prs,cellid,varargin{:})
g = prs.Results;

% Load event structure
event_type = lower(g.event_type(1:4));
switch event_type
    case 'stim'
        
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');
        catch ME
            disp('There was no stim protocol for ths session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');
        catch ME
            disp('There was no behavioral protocol for ths session.')
            error(ME.message)
        end
        
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

% Calculate time stamps
tno = length(g.abstimes);
reltimes = nan(1,tno);
if isequal(event_type,'tria') && ~isequal(event,'TrialStart')
    for k = 1:tno
        pr = g.abstimes(k) - (VE.(g.event) + VE.TrialStart);     % trial events other than TrialStart are stored relative to TrialStart
        pr = pr(pr>0);
        reltimes(k) = pr(end);
    end
else
    for k = 1:tno
        pr = g.abstimes(k) - VE.(g.event);
        pr = pr(pr>0);
        reltimes(k) = pr(end);
    end
end