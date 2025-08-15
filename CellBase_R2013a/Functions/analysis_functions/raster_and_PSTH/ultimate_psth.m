function stats = ultimate_psth(cellid,spt,varargin)
%ULTIMATE_PSTH   Peri-stimulus time histogram.
%   [PSTH SPSTH SPSTH_SE] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates peri-stimulus time histogram (PSTH) for the cell passed in
%   CELLID. Smoothed PSTH (SPSTH) and SE of smoothing (SPSTH_SE) are also
%   returned.
%
%   [PSTH SPSTH SPSTH_SE TAGS] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   returns partition tags (TAGS) corrsponding to PSTHs when trials are 
%   partitioned; see PARTITION_TRIALS.
%
%   [PSTH SPSTH SPSTH_SE TAGS SPT] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   returns the bin raster (SPT); see STIMES2BINRASTER.
%
%   [PSTH SPSTH SPSTH_SE SPT TAGS STATS] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates and returns test results for significant firing rate changes
%   after the event (see PSTH_STATS for details).
%
%   ULTIMATE_PSTH is also capable of using two different events for the
%   periods before and after 0, usefull for statistical testing with a
%   baseline period aligned to a different event than the test period (see
%   below and PSTH_STATS).
%
%   Mandatory input arguments:
%       CELLID: defines the cell (see CellBase documentation)
%       EVENT: the event to which the PSTH is aligned; if EVENT is a cell
%           array of two strings, the first event is used for the PSTH 
%           and binraster before 0 and the second event is used for the 
%           PSTH and binraster after 0; if EVENT is a function handle, the
%           function is called for CELLID to define the aligning event
%           (dynamic event definition)
%       EVENT_TYPE: the type of event, 'stim' or 'trial'
%       WINDOW: window for calculation relative to the event in seconds
%
%   Default behavior of ULTIMATE_PSTH can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'dt', 0.001 - time resolution in seconds
%       'sigma', 0.02 - smoothing kernel for the smoothed PSTH, in seconds
%       'margin',[-0.01 0.01] margins for PSTH calculation to get rid of 
%           edge effect due to smoothing
%       'event_filter', 'none' - filter light-stimulation trials; see
%           FILTERTRIALS for implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%       'maxtrialno', 5000 - maximal number of trials included; if ther are
%           more valid trials, they are randomly down-sampled
%       'parts', 'all' - partitioning the set of trials; input to
%           PARTITION_TRIALS, see details therein (default, no
%           partitioning)
%       'isadaptive', 1 - 0, classic PSTH algorithm is applied; 1, adaptive
%           PSTH is calculated (see APSTH); 2, 'doubly adaptive' PSTH
%           algorithm is used (see DAPSTH)
%   	'baselinewin', [-0.25 0] - limits of baseline window for
%           statistical testing (see PSTH_STATS), time relative to 0 in 
%           seconds
%   	'testwin', [0 0.1] - limits of test window for statistical testing
%           (see PSTH_STATS), time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals in PSTH_STATS; in
%           proportion of the peak-baseline difference (see PSTH_STATS)
%       'display', false - controls plotting
%
%   See also PSTH_STATS, STIMES2BINRASTER, BINRASTER2PSTH, BINRASTER2APSTH,
%   APSTH, VIEWCELL2B, PARTITION_TRIALS and FILTERTRIALS.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   07-May-2012

%   Edit log: BH 7/5/12, 8/12/12, 8/27/12

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@(s)ischar(s)|...
    (iscellstr(s)&isequal(length(s),2))|...
    isa(s,'function_handle'))   % reference event
addParameter(prs,'window',[-0.1 0.1],...
    @(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'event_filter','none',@(s)ischar(s)|iscellstr(s))   % filter events based on properties
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParameter(prs,'sigma',0,@isnumeric)     % smoothing kernel for the smoothed PSTH
% addParamValue(prs,'margin',[-0.1 0.1])  % [-0.1,0.1] margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'parts','all')   % partition trials
addParameter(prs,'isadaptive',true,@(s)islogical(s)|ismember(s,[0 1 2]))   % use adaptive PSTH algorithm
addParameter(prs,'baselinewin',[-0.25 0],@(s)isnumeric(s)&isequal(length(s),2))  % [-0.25 0],time window relative to the event for stat. testing, in seconds
addParameter(prs,'testwin',[0 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % [0 0.1] time window relative to the event for stat. testing, in seconds
addParameter(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>0&s<1)   % threshold used to assess interval limits in PSTH_STATS
addParameter(prs,'display',true,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;
if nargout > 2   % statistics will only be calculted if asked for
    g.dostats = true;
else
    g.dostats = false;
end

% Load event structure
switch g.event_type
    case 'stim'
        
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');   % load events
            VS = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes
        catch ME
            disp('There was no stim protocol for this session.')
            error(ME.message)
        end
        
    case 'tria'
        
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');   % load events
            VS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
        catch ME
            disp('There was no behavioral protocol for this session.')
            error(ME.message)
        end
        
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

event=g.event;
% Events, time, valid trials
margin = g.sigma * 3;     % add an extra margin to the windows
time = g.window(1)-margin:g.dt:g.window(2)+margin;   % time base array
if iscell(event)   % different ref. event for tst and baseline window
    event1 = event{1};
    event_pos1 = findcellstr(VS.events(:,1),event1);
    event2 = event{2};
    event_pos2 = findcellstr(VS.events(:,1),event2);
    if event_pos1 * event_pos2 == 0
        error('Event name not found');
    end
    stimes1 = VS.event_stimes{event_pos1};
    stimes2 = VS.event_stimes{event_pos2};
    valid_trials1 = filterTrials(cellid,'event_type',g.event_type,'event',event1,...
        'event_filter',g.event_filter,'filterinput',g.filterinput);
    valid_trials2 = filterTrials(cellid,'event_type',event_type,'event',event2,...
        'event_filter',g.event_filter,'filterinput',g.filterinput);
    if ~isequal(valid_trials1,valid_trials2)
        error('Valid trials should be the same for baseline and test period.')
    else
        valid_trials = valid_trials1;
    end
else
    if isa(event,'function_handle')
        event = feval(event,cellid);   % dynamic event definition
    end
    event_pos = findcellstr(VS.events(:,1),event);
    if event_pos == 0
        error('Event name not found');
    end
    stimes = VS.event_stimes{event_pos};
    valid_trials = filterTrials(cellid,'event_type',g.event_type,'event',g.event,...
        'event_filter',g.event_filter,'filterinput',g.filterinput);
end

% Calculate bin rasters
% if iscell(event)   % different ref. event for tst and baseline window
%     spt1 = stimes2binraster(stimes1(valid_trials),time,g.dt);
%     spt2 = stimes2binraster(stimes2(valid_trials),time,g.dt);
%     spt = [spt1(:,time<0) spt2(:,time>=0)];
% else
%     spt = stimes2binraster(stimes(valid_trials),time,g.dt);
% end

% Partition trials
[COMPTRIALS, ~] = partition_trials(VE,g.parts);

% PSTH
switch g.isadaptive
    case {0,false}
        [~, spsth, ~] = binraster2psth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
    case {1, true}
        [~, spsth, ~] = binraster2apsth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
    case 2
        [~, spsth, ~] = binraster2dapsth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
end


% Output statistics
if g.dostats
    stats = psth_stats(cellid,spt,spsth,'dt',g.dt,'win',g.window,...
        'sigma',g.sigma,'baselinewin',[g.window(1) 0],'testwin',[0.001 g.window(2)],'display',0,...
        'relative_threshold',g.relative_threshold);
end