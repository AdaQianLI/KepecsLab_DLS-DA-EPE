function [Eff,L,J] = reliability_latency_jitter(cellid,stats,varargin)

%RELIABILITY_LATENCY_JITTER   Reliability, latency and jitter of spiking.
%   [R L J] = RELIABILITY_LATENCY_JITTER(CELLID) calculates the latency of
%   spiking of a cell (CELLID) after an event (default, 'PulseOn') (L). It
%   also returns the standard deviation of spike times after finding the
%   period of evoked spikes using ULTIMATE_PSTH (J). Reliability of evoking
%   spikes (proportion of events followed by evoked spikes) is returned in
%   R.
%
%   [R L J B M] = RELIABILITY_LATENCY_JITTER(CELLID,...) also returns
%   baseline firing rate (B) and maximal firing rate in the test window
%   (M).
%
%   [R L J B M A1 A2] = RELIABILITY_LATENCY_JITTER(CELLID,...) returns the
%   start and end point of detected stimulation period (see
%   ULTIMATE_PSTH and FINDSTIMPERIOD).
%
%   [R L J B M A1 A2 SPNO] = RELIABILITY_LATENCY_JITTER(CELLID,...) returns
%   the number of spikes detected for each valid trial (see EFFICIENCY).
%   SPNO: number of spikes for all valid trials
%   [R L J B M A1 A2 SPNO H] = RELIABILITY_LATENCY_JITTER(CELLID,...)
%   returns figure handles for PSTH and raster plot in the structure H (if
%   'display' is set to true).
%
%   Optional parameter-value pairs, with default values:
%       'event_type', 'stim' - type of the aligning event, 'stim' (stimulus
%           events) or 'trial' (trial events)
%       'event', 'PulseOn' - aligning event
%       'event_filter', 'none' - filter trials; see FILTERTRIALS for 
%           implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%       'window', [-0.005 0.01] - time window relative to the event, in
%           seconds; for determining the range of the PSTH, see
%           ULTIMATE_PSTH
%       'isadaptive', 1 - 0, classic PSTH algorithm is applied; 1, adaptive
%           PSTH is calculated (see APSTH); 2, 'doubly adaptive' PSTH
%           algorithm is used (see DAPSTH)
%   	'baselinewin', [-0.005 0] - limits of the baseline window for PSTH 
%           statistics (see PSTH_STATS), time relative to 0 in seconds
%   	'testwin', [0 0.01] - limits of the test window for PSTH statistics
%           (see PSTH_STATS), time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals in PSTH_STATS; in
%           proportion of the peak-baseline difference (see PSTH_STATS); it
%           determines the window within which spikes will be considered
%           'evoked'
%       'jitterdefinition', 'all' - control how jitter is defined: 
%           'all' - SD of spike times of all spikes
%           'burst' - only first spikes are included in each trial
%       'display', false - controls plotting
%
%   Examples:
%   [E L J] = reliability_latency_jitter(cellid,'event','BurstOn');
%
%   fi = struct('BurstNPulse',20);
%   [E I J] = reliability_latency_jitter(cellid,...
%       'event_filter','BurstNPulse_maxPower','filterinput',fi);
%
%   [E_burston L_burston J_burston B_burston M_burston Astart Aend] = ...
%        reliability_latency_jitter(cellid,'event','PulseOn');
%
%   See also ULTIMATE_PSTH, EXTRACTSEGSPIKES and EFFICIENCY.
%% Wpa' - p value for Mann-Whitney test for significant activation
%% 'Wpi' - p value for Mann-Whitney test for significant inhibition

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   15-June-2013

%   Edit log: BH 7/15/13

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'event_type','stim',...
    @(s)ischar(s)&&(strncmp(s,'stim',4)||strncmp(s,'trial',4)))   % event type - stim or trial
addParameter(prs,'event','PulseOn_chr2',@ischar)   % reference event
addParameter(prs,'event_filter','none',@ischar)   % event filter
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'window',[-0.1 0.1],...
    @(s)isnumeric(s)&isequal(length(s),2))  % [-0.005 0.01],time window relative to the event, in seconds
addParameter(prs,'dt',0.001,@isnumeric)   
addParameter(prs,'sigma',0.001,@isnumeric)  % time resolution of the binraster, in seconds
addParameter(prs,'isadaptive',1,@(s)islogical(s)|ismember(s,[0 1 2]))   % use adaptive PSTH algorithm
addParameter(prs,'baselinewin',[-0.1 0],@(s)isnumeric(s)&isequal(length(s),2))  % [-0.005 0],time window relative to the event for stat. testing, in seconds
addParameter(prs,'testwin',[0 0.01],@(s)isnumeric(s)&isequal(length(s),2))  % [0 0.01] time window relative to the event for stat. testing, in seconds
addParameter(prs,'relative_threshold',0.2,...
    @(s)isnumeric(s)&isequal(length(s),1)&s>0&s<1)  % 0.5, relative threshold for peak detection, see ULTIMATE_PSTH
addParameter(prs,'jitterdefinition','all',...
    @(s)ischar(s)|ismember(s,{'all','burst'}))   % controls the definition of 'jitter'
addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;
% H = struct('H_psth',{},'H_raster',{});   % initialize structure for figure handles



% Filter events
valid_trials = filterTrials(cellid,'event_type',g.event_type,'event',g.event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);


lim_act1 = stats.activation_start;  % window start for evoked spikes
lim_act2 = stats.activation_end*2;   % window end for evoked spikes


Base = stats.baseline;  % baseline FR
Max = stats.maxvalue;  % peak FR
Min=stats.minvalue;
lim_inh1=stats.inhibition_start;
lim_inh2=stats.inhibition_end;


if (Max-Base)>(Base-Min)
L = stats.activation_peak;   % LATENCY
tsegs_evoked = rel2abstimes(cellid,[lim_act1 lim_act2],'event_type',g.event_type,'event',g.event,'valid_trials',valid_trials);   % convert times rel. to the event to absolute spike times
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find evoked spikes
else
L=stats.inhibition_peak;
tsegs_evoked = rel2abstimes(cellid,[lim_inh1 lim_inh2],'event_type',g.event_type,'event',g.event,'valid_trials',valid_trials);   %  convert times rel. to the event to absolute spike times
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find evoked spikes
end
% Reliability
Eff = efficiency_spike(cellid,selts_evoked,'event_type',g.event_type,'event',g.event,'valid_trials',valid_trials);   % efficiency of evoking spikes

% Standard deviation of evoked spikes (jitter)
switch g.jitterdefinition
    case 'all'
        relevokedtimes = abs2reltimes(cellid,selts_evoked,'event_type',g.event_type,'event',g.event);  % convert absolute spike times to times rel. to the event
    case 'burst'
        relevokedtimes = abs2reltimes(cellid,selts_evoked1st,'event_type',g.event_type,'event',g.event);  % convert absolute spike times to times rel. to the event
end
J = std(relevokedtimes);   % JITTER

end