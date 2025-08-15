function H = PlotStimPeriod(cellid,varargin)
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
% addRequired(prs,'reltimes',@(s)isnumeric(s)|iscell(s))  % time stamps to convert
addParameter(prs,'event_type','stim',...
    @(s)ischar(s)&&(strncmp(s,'stim',4)||strncmp(s,'trial',4)))   % event type ('stim' or 'trial')
addParameter(prs,'event','PulseOn',@ischar)   % reference event
addParameter(prs,'valid_trials','all')   % valid trials - use all trials by default
addParameter(prs,'window',[-0.02 0.2],...  % [-0.005 0.04]
    @(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'limwindow',[],...  % [-0.005 0.04]
    @(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'dt',0.001,@isnumeric)   % 0.0005, time resolution of the binraster, in seconds
addParameter(prs,'margin',[-0.01 0.01])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;

if isnan(g.limwindow(1)) || isnan(g.limwindow(2))
    g.limwindow(1) = 0.001;   % no activation detected
    g.limwindow(2) = 0.006;
end

% % % Directories
% % global DATAPATH
% % fs=filesep;
% % DATAPATH = uigetdir;
% % resdir = [DATAPATH 'Taggedprop'];

% Set parameters and load CellBase variables
datapath = getpref('cellbase','datapath');
[subjectID, sessionID] = cellid2tags(cellid);
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
if g.display
    SEvent = 'BurstOff';
    FNum = 2;
    parts = 'all';
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent = {{'BurstOff'}};
    ShEvColors = hsv(length(ShEvent{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    
    % Plot raster plot and PSTH for 'PulseOn'
    set(gcf,'renderer','painters')   % temporarily change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',g.event,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',g.window,'dt',g.dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on')
    pause(0.05)   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

% PSTH
dtt = g.dt * 1000;   % resolution of bin raster in ms
wn = g.window * 1000;   % window boundaries in ms
margin = g.margin * 1000;   % margin in ms

H = newpsth(spt,spt,dtt,wn,margin,g.limwindow);


%% save H handle plot
cellidt=cellid;
cellidt(cellidt=='.') = '_';
resdir = [datapath '/' subjectID '/' sessionID 'Summary/'];
if ~isdir(resdir)
    mkdir(resdir)
end

        fnm = [resdir cellidt '_' 'SDF.fig'];   % save
        saveas(H,fnm);
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
function H= newpsth(spt_baseline,spt_test,dt,win,margin,limwindow)

% Trial number and epoch length
[tno_baseline,tl] = size(spt_baseline);
[tno_test,tl] = size(spt_test);

% Pre-stimulus time window to consider for null hypothesis
stm = abs(win(1)+margin(1)) / dt;

% Merged spike train
sptb = spt_baseline(:,1:stm);
[x0,allspks_baseline] = find(sptb);
ts_baseline = sort(allspks_baseline)';

sptt = spt_test(:,stm+1:end);
[x0,allspks_test] = find(sptt);
ts_test = stm + sort(allspks_test)';
ts = [ts_baseline ts_test];

% Calculate adaptive SDF with variable Gaussian Kernel
prob = [sum(sptb)/tno_baseline sum(sptt)/tno_test]  / dt;  % prob. if 1 ms bins; dt is measured in ms!
spno = length(ts);
agvd = zeros(1,tl);
for t = 1:spno
    spi = ts(t);
    tspt = zeros(1,tl);
    tspt(spi) = 1;
    wbh = gausswin(9,prob(spi)*50);   % kernel
    wbh = wbh / sum(wbh);
    if spi <= stm
        agvd = agvd + filtfilt(wbh,1,tspt) / tno_baseline;   % convolution from both directions
    else
        agvd = agvd + filtfilt(wbh,1,tspt) / tno_test;
    end
end
ppsth_aconv = agvd / dt * 1000;   % SDF
psth_aconv = ppsth_aconv((stm+1+win(1)/dt):(stm+1+win(2)/dt));
% psth_aconv_ave=normr(psth_aconv);
% prob = prob((stm+1+win(1)/dt):(stm+1+win(2)/dt));
% stn = abs(win(1)) / dt;

% Plot SDF
time = win(1):dt:win(end);
lim1=round(limwindow(1)/dt*1000);  %change s to ms
lim2=round(limwindow(2)/dt*1000);  

   win_low  = find(time <= lim1,1,'last');
   win_high = find(time >= lim2,1,'first');
    H = figure;
    subplot(3,1,1);
    plot(time,psth_aconv,'k');
    hold on;
    plot(time(win_low:win_high),psth_aconv(win_low:win_high),'Color',[0.8 0 0],'LineWidth',2);
    xlim([time(1) time(end)]);
    xlabel('time [ms]')
    ylabel('firing rate');
    subplot(3,1,2);
    cum_psth = cumsum(psth_aconv/sum(psth_aconv));
    plot(time,cum_psth);
    xlim([time(1) time(end)]);
    ylim([0,1]);
    subplot(3,1,3);
    cdfplot(psth_aconv);


% Activation/inhibition time
% baseline_prob = mean(prob(1:stn)) * 1000;  % spikes/sec (was spikes/bin before)
% nst = length(psth_aconv);   % changing this allows further restricting the window
% maxafter = max(psth_aconv(stn+1:nst));
% minafter=min(psth_aconv(stn+1:nst));
% if (maxafter - baseline_prob) < (baseline_prob-minafter)   % putative activation, if firing goes above baseline
%     mininx = stn + find(psth_aconv(stn+1:nst)==minafter,1,'first');   % minimun firing
%     thr = baseline_prob - (baseline_prob - minafter)/2;   % 0.5 default
%     pas = valuecrossing(time(stn+1:mininx),psth_aconv(stn+1:mininx),thr,'down');
%     pas_inx = valuecrossing(stn+1:mininx,psth_aconv(stn+1:mininx),thr,'down');
%     if  isempty(pas)
%         pas = time(stn+1);
%         pas_inx = stn + 1;
%     end
%     pas_inx = round(pas_inx(end));
%     activation_start = pas(end);
%     pae = valuecrossing(time(mininx:nst),psth_aconv(mininx:nst),thr,'up');
%     pae_inx = valuecrossing(mininx:nst,psth_aconv(mininx:nst),thr,'up');
%      if isempty(pae)
%         pae = time(nst);
%         pae_inx = nst;
%      end
%     pae_inx = round(pae_inx(1));
%     activation_end = pae(1);
%     activation_time = activation_end - activation_start;
%     activation_peak = time(mininx) - time(stn+1);   
% % %     activation_start = NaN;
% % %     activation_end = NaN;
% % %     activation_peak = NaN;
% % %     activation_time = 0;   % if firing does not go above baseline
% % %     pas_inx=1;
% % %     pae_inx=length(time);
% else
%     maxinx = stn + find(psth_aconv(stn+1:nst)==maxafter,1,'first');   % maximal firing
%     thr = baseline_prob + (maxafter - baseline_prob) *0.2;   % 0.5 default
%     pas = valuecrossing(time(stn+1:maxinx),psth_aconv(stn+1:maxinx),thr,'up');
%     pas_inx = valuecrossing(stn+1:maxinx,psth_aconv(stn+1:maxinx),thr,'up');
%     if isempty(pas)
%         pas = time(stn+1);
%         pas_inx = stn + 1;
%     end
%     pas_inx = round(pas_inx(end));
%     activation_start = pas(end);   % last crossing of one and a half-baseline probability before maximum
%     pae = valuecrossing(time(maxinx:nst),psth_aconv(maxinx:nst),thr,'down');
%     pae_inx = valuecrossing(maxinx:nst,psth_aconv(maxinx:nst),thr,'down');
%     if isempty(pae)
%         pae = time(nst);
%         pae_inx = nst;
%     end
%     pae_inx = round(pae_inx(1));
%     activation_end = pae(1);   % first crossing of one and a half-baseline probability after maximum
%     activation_time = activation_end - activation_start;
%     activation_peak = time(maxinx) - time(stn+1);    % peak time of activation
% end

% % Plot
%     figure(H)
%     subplot(3,1,1); 
%     hold on;
%     try
%     H=
%     end
%     y_lim = ylim;
%     ylim([0 y_lim(2)])
%     xlabel('time [ms]')
%     ylabel('firing rate');
%     setappdata(H,'activation_start',activation_start)   % store the variables with the figure
%     setappdata(H,'activation_end',activation_end)
%     setappdata(H,'activation_peak',activation_peak)
%     setappdata(H,'activation_time',activation_time)
%     setappdata(H,'baseline',baseline_prob)
%     setappdata(H,'maxafter',maxafter)

end
