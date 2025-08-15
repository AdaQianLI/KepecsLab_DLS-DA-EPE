function H=findStimPeriod(cellid,varargin)
%FINDSTIMPERIOD   Find time period of stimulated spikes.
%   [I1 I2] = FINDSTIMPERIOD(CELLID) seeks for an increase of firing rate
%   after an event. First, the program calculates Spike Density Function by
%   convolving the event-aligned raster plots with a variable Gaussian
%   window (see SOMPSTH_CALL2). Second, it finds maximal firing as maximum
%   SDF within a window from the event. Baseline firing is determined by
%   mean pre-event firing probability. Next, the time course of activation
%   is assessed by half-baseline crossings before and after the maximum.
%   This temporal window of activation is then returned in I1 (start) and
%   I2 (end timestamp).
%
%   [I1 I2 P] = FINDSTIMPERIOD(CELLID) returns the peak time of activation
%   (in seconds) as third output argument (P).
%
%   [I1 I2 P H] = FINDSTIMPERIOD(CELLID) plots the adaptive SDF and returns
%   the figure handle in H. Activation parameters are stored with the
%   figure as application data.
%
%   FINDSTIMPERIOD accepts the following parameter, value pairs as optional
%   input arguments (with default values):
%       'event', 'PulseOn' - the event to which the window is locked
%   	'window', [-0.005 0.01] - extent of baseline and test period
%           relative to the event in seconds; in seconds
%   	'margins', [-0.01 0.01] - margin in seconds used for extending the
%           'window' for adaptive SDF calculation to reduce convolution edge
%           effects; the SDF is later restricted to the 'window'
%       'dt', 0.0005 - time resolution of the bin raster; in seconds
%       'display', false - control of plotting event-locked raster plot
%       'valid_trials', 'all' - filter for trials; default: include all
%           trials; logical array or index set
%
%   See also ULTIMATE_PSTH, FINDSEGS3 and ABS2RELTIMES.

%   Edit log: BH 4/25/12, 5/4/12

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'event_type','stim',...
    @(s)ischar(s)&&(strncmp(s,'stim',4)||strncmp(s,'trial',4)))   % event type ('stim' or 'trial')
addParameter(prs,'event','PulseOn',@ischar)   % reference event
addParameter(prs,'valid_trials','all')   % valid trials - use all trials by default
addParameter(prs,'window',[-0.02 0.2],...  % [-0.005 0.04]
    @(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'dt',0.001,@isnumeric)   % 0.0005, time resolution of the binraster, in seconds
addParameter(prs,'margin',[-0.01 0.01])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;

% % % Directories
% % global DATAPATH
% % fs=filesep;
% % DATAPATH = uigetdir;
% % resdir = [DATAPATH 'Taggedprop'];

% Set parameters and load CellBase variables
% datapath = getpref('cellbase','datapath');
% [subjectID, sessionID] = cellid2tags(cellid);
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
SE = loadcb(cellid,'StimEvents');
epoch_pos = findcellstr(ST.events(:,1),g.event);
if epoch_pos == 0
    error('Event not found.');
end
stimes = ST.event_stimes{epoch_pos};
time = (g.window(1)+g.margin(1)):g.dt:(g.window(2)+g.margin(2));

% Valid trials
valid_trials = parseValidTrials(SE,g.event,g.valid_trials);

% Calculate bin raster
spt = stimes2binraster(stimes(valid_trials),time,g.dt);

% Set input arguments for raster plot and PSTH
% if g.display
%     SEvent = 'BurstOff';
%     FNum = 2;
%     parts = 'all';
%     sigma = 0.001;
%     PSTHstd = 'on';
%     ShEvent = {{'BurstOff'}};
%     ShEvColors = hsv(length(ShEvent{1}));
%     ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
%     
%     % Plot raster plot and PSTH for 'PulseOn'
%     set(gcf,'renderer','painters')   % temporarily change renderer because OpenGL locks the plot which result an error in legend layout handling
%     viewcell2b(cellid,'TriggerName',g.event,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
%         'FigureNum',FNum,'eventtype','stim','window',g.window,'dt',g.dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
%         'EventMarkerWidth',0,'PlotZeroLine','on')
%     pause(0.05)   % if reset the renderer two early, the same error occurs
%     set(gcf,'renderer','opengl')   % reset renderer
% end

% PSTH
dtt = g.dt * 1000;   % resolution of bin raster in ms
wn = g.window * 1000;   % window boundaries in ms
margin = g.margin * 1000;   % margin in ms
H = newpsth(spt,dtt,wn,margin);
% activation_start = activation_start / 1000;   % convert back to seconds
% activation_end = activation_end / 1000;
% activation_peak = activation_peak / 1000;
% activation_time = activation_time / 1000;

% %% save H handle plot
% cellidt=cellid;
% cellidt(cellidt=='.') = '_';
% resdir = [datapath '/' subjectID '/' sessionID 'PSTH/'];
%         fnm = [resdir cellidt '_' 'SDF.fig'];   % save
%         saveas(H,fnm);
end
% % if nargout > 6
% % %%%%%%calculate SALT
% % sum_evokedspike=sum(evokedspike);
% % peak_idx=find(sum_evokedspike==max(sum_evokedspike));
% % spt_baseline=spt;
% % spt_baseline(:,91:end)=[];
% % spt_baseline(:,31:40)=[];
% % spt_baseline(:,1:10)=[];
% % [row_spt col_spt]=size(spt_baseline);
% % spt_baseline=reshape(spt_baseline,row_spt/5,col_spt*5);
% % spt_evoked=evokedspike(1:100,1:10);
% % [row_evoked col_evoked]=size(spt_evoked);
% % window=col_evoked*g.dt;
% % [p I] = salt2(spt_baseline,spt_evoked,g.dt,window);
% % end
% -------------------------------------------------------------------------

