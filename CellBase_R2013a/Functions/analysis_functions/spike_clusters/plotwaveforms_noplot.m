function [R,p_v_width,wv] = plotwaveforms_noplot(cellid,stats,ADBitVolts,varargin)
%PLOTWAVEFORMS   Plot spontaneous and/or light-evoked spike shape.
%   PLOTWAVEFORMS(CELLID) plots the waveform of spontaneous and putative
%   light-evoked spikes. It finds the latency of putative light-evoked
%   spikes as peaks of firing rate after light pulses. The period for
%   putative light-evoked spikes is defined by half-hight crossings around
%   the peak (see FINDSTIMPERIOD for details). If there is no detectable
%   activation, 1 and 6 ms are used, which include typical latencies of
%   light-evoked spikes. Spontaneous spikes are extracted from 2 s periods
%   before light bursts. Only 2000 randomly chosen spikes are plotted. This
%   number can be changed using the 'MAXNUM', VALUE optional input
%   parameter-value pair.
%
%   Two optional input parameter value pairs control whether to plot both
%   spontaneous and evoked waveforms or only one of those:
%   HS = PLOTWAVEFORMS(CELLID,'SPONT',TRUE,'EVOKED',TRUE) plots both
%   waveforms and a third plot to compare them. Setting either of the input
%   values to 'false' will prevent plotting that type. Handles of the
%   resulting figures are returned in the output structure HS with the
%   following fieldnames: HS.H_spont for spontaneous waveforms, HS.H_evoked
%   for light-evoked waveforms and HS.H_compare for the compound figure.
%
%   HS = PLOTWAVEFORMS(CELLID,'CORRELATION','TRUE') also returns waveform
%   correlation for the largest channel (see SPIKESHAPECORR) in HS.R field
%   of the output structure.
%
%   Optional input parameter-value paris with default values:
%       'spont', true - plot spontaneous waveform
%       'evoked', true - plot light-evoked waveform
%       'maxnum', 2000 - maximum number of spikes to plot
%       'correlation', false - return waveform correlation
%       'stim_period', [] - start and end of stimulation period after each
%           light pulse
%
%   Output (fields of HS struct):
%       H_spont - figure handle for spontaneous spike waveform
%       H_evoked - figure handle for evoked spike waveform
%       H_compare - figure handle for comparison plot
%       R - spike shape correlation of spont. and evoked waveform
%       activation_start - start of the detected stimulation period
%       activation_end - end of the detected stimulation period
%
%   See also FINDSTIMPERIOD, ABS2RELTIMES, EXTRACTSEGSPIKES,
%   EXTRACTSPIKEWAVEFORMS and SPIKESHAPECORR.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   07-May-2012

%   Edit log: BH 5/7/12

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'event','PulseOn_chr2',@ischar)   % reference event
addParameter(prs,'event_filter','none',@ischar)   % event filter
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'event_type','stim',...
    @(s)ischar(s)&&(strncmp(s,'stim',4)||strncmp(s,'trial',4)))   % event type - stim or trial
addParameter(prs,'spont',true,@(s)islogical(s)|ismember(s,[0 1])) % plot spontaneous waveform
addParameter(prs,'evoked',true,@(s)islogical(s)|ismember(s,[0 1])) % plot light-evoked waveform
addParameter(prs,'maxnum',2000,@isnumeric)   % maximal number of spikes to plot
addParameter(prs,'correlation',true,@(s)islogical(s)|ismember(s,[0 1])) % calculate waveform correlation
addParameter(prs,'stim_period',[],@isnumeric)   % start and end point for the interval of light-evoked spikes
addParameter(prs,'spon_window',[-0.1 0.1],...
    @(s)isnumeric(s)&isequal(length(s),2)) 
addParameter(prs,'display',true,@(s)islogical(s)|ismember(s,[0 1])) % plot spontaneous waveform

parse(prs,cellid,varargin{:})
g = prs.Results;

% % global DATAPATH
% % DATAPATH = uigetdir;
% % resdir = [DATAPATH 'Taggedprop'];


% Evoked spikes

%     ST = loadcb(cellid,'STIMSPIKES');  % load stimulus spikes (prealigned)
%     if isequal(findcellstr(ST.events(:,1),g.event),0)  % if prealignment to stim. period has not been saved yet
%         if isempty(g.stim_period)   % latency of stimulated spikes
%             [lim1,lim2,~,~,~,~,H] = findStimPeriod(cellid);   % find putative stimulated period
%             if isnan(lim1) || isnan(lim2)
%                 lim1 = 0.001;   % no activation detected
%                 lim2 = 0.006;
%             end
%         else  % stimulation period passed as input
%             lim1 = g.stim_period(1);
%             lim2 = g.stim_period(2);
%         end
%         tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'stim',g.event);   % convert period to epochs relative to pulses
%         selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find putative stimualated spikes
%     else  % if prealignment to stim. period has already been saved
%         trigger_pos = findcellstr(ST.events(:,1),g.event);
%         pse = ST.event_stimes{trigger_pos};
%         selts_evoked = rel2abstimes(cellid,pse,'event_type','stim','event',g.event);  % convert to absolute spike times
%         lim1 = ST.events{trigger_pos,4}(1);  % boundaries of stim. period used for prealigning spikes
%         lim2 = ST.events{trigger_pos,4}(2);
%     end

% Filter events
valid_trials = filterTrials(cellid,'event_type',g.event_type,'event',g.event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);


% lim_act1 = stats.activation_start;  % window start for evoked spikes
% lim_act2 = stats.activation_end;   % window end for evoked spikes


lim_act1 = 0;  % window start for evoked spikes
lim_act2 = 0.02;   % window end for evoked spikes

Base = stats.baseline;  % baseline FR
Max = stats.maxvalue;  % peak FR
Min=stats.minvalue;
lim_inh1=stats.inhibition_start;
lim_inh2=stats.inhibition_end;


if (Max-Base)>(Base-Min)
    stimwin=[lim_act1 lim_act2];
tsegs_evoked = rel2abstimes(cellid,stimwin,'event_type',g.event_type,'event',g.event,'valid_trials',valid_trials);   % convert times rel. to the event to absolute spike times
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find evoked spikes
selts_evoked=unique(selts_evoked);
else
stimwin=[lim_inh1 lim_inh2];
tsegs_evoked = rel2abstimes(cellid,stimwin,'event_type',g.event_type,'event',g.event,'valid_trials',valid_trials);   %  convert times rel. to the event to absolute spike times
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find evoked spikes
selts_evoked=unique(selts_evoked);
end


    wave_evoked = extractSpikeWaveforms(cellid,selts_evoked,'chans','all');  % get waveforms for the extracted spikes
    wave_evoked = wave_evoked* ADBitVolts * 1e6;
    
    % Downsample
    nm_evoked = size(wave_evoked,1);
    rp = randperm(min(g.maxnum,nm_evoked));
    weds = wave_evoked(rp,:,:);   % downsample for plotting
    
    % Average waveforms
    mean_evoked = squeeze(nanmean(wave_evoked,1));


% Spontaneous spikes
if g.spont || g.correlation
    tsegs_spont = rel2abstimes(cellid,g.spon_window,'event_type','stim','event',g.event);   % extract 2s periods before bursts
    selts_spont = extractSegSpikes(cellid,tsegs_spont);     % extract spontaneous spikes
    selts_spont=unique(selts_spont);
    
    wave_spont = extractSpikeWaveforms(cellid,selts_spont,'chans','all');    % get waveforms for the extracted spikes
    wave_spont = wave_spont* ADBitVolts * 1e6;
    %     wave_prop=calcwaveformfeatures(wave_spont);
    % Downsample
    nm_spont = size(wave_spont,1);
    rp = randperm(min(g.maxnum,nm_spont));
    wsds = wave_spont(rp,:,:);
    
    % Average waveforms
    mean_spont = squeeze(nanmean(wave_spont,1));
end

% % Plot light-evoked waveforms
% if g.evoked
%     out.H_evoked = figure('Position',[624 126 1092 852]);
%     H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
%     for sp = 1:4
%         hold(H(sp),'on')
%         plot(H(sp),transpose(squeeze(weds(:,sp,:))))
%     end
%     title('Light-evoked spike shape')
% end
% % % filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_Light_evoked_spike_shape.fig'];
% % % saveas(gcf, filename);
% % Plot spontaneous waveforms
% if g.spont
%     out.H_spont = figure('Position',[624 126 1092 852]);
%     H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
%     for sp = 1:4
%         hold(H(sp),'on')
%         plot(H(sp),transpose(squeeze(wsds(:,sp,:))))
%     end
%     title('Spontaneous spike shape')
% end
% % filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_Spontaneous_spike_shape.fig'];
% % saveas(gcf, filename);

% Compare waveforms
if g.evoked && g.spont
%     out.H_compare = figure('Position',[624 126 1092 852]);
%     H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
if g.display
    if (~isnan(mean_evoked(1,1)) && ~isnan(mean_spont(1,1)))
        for sp = 1:4
        subplot(2,2,sp);        
        hold on;
        H{sp,1}=plot(transpose(squeeze(wsds(:,sp,:))),'Color',[0.9 0.9 0.9]);
        H{sp,2}=plot(transpose(squeeze(weds(:,sp,:))),'Color',[0 150 255]/255);
        H{sp,3}=plot(transpose(mean_spont(sp,:)),'Color','k','LineWidth',2);  
        H{sp,4}=plot(transpose(mean_evoked(sp,:)),'Color',[0 50 255]/255,'LineWidth',2);
        xlim([1 size(wsds,3)]);
        end
    else
        for sp = 1:4
        subplot(2,2,sp); 
        hold on
        H{sp,1}=plot(ones(1,32),'.w');
        H{sp,2}=plot(ones(1,32),'.w');
        H{sp,3}=plot(ones(1,32),'.w');
        H{sp,4}=plot(ones(1,32),'.w');
        end
    end
    title('Compare spont. and light-evoked spike shape')
% %     filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_Light_Spon_spike_shape.fig'];
% %     saveas(gcf, filename);
end

        mx = maxchannel(wave_spont);     % mx: largest channel
        if size(squeeze(wave_spont(:,mx,:)),1)>1
        wv = nanmean(squeeze(wave_spont(:,mx,:)));
        else
        wv = squeeze(wave_spont(:,mx,:));    
        end   
% Spike shape correlation
try
        mx = maxchannel(wave_spont);     % mx: largest channel
        mnmx_spont = nanmean(squeeze(wave_spont(:,mx,:)));
        mnmx_evoked = nanmean(squeeze(wave_evoked(:,mx,:)));
        sr = 32552;     % DigiLynx sampling rate
        rng = round(0.00075*sr);    % number of data points in 750 us (default censored period of DigiLynx)
        pr = corrcoef(mnmx_spont(1:rng),mnmx_evoked(1:rng));
        R = pr(1,2);
catch
        R = nan;
end
if g.display
Y=text(20,8,['R = ',num2str(R)],'Color',[1 0 0],'FontSize',8);
end
end

if ~isempty(mnmx_spont)
 [~,idx1]=max(mnmx_spont);
 [~,idx2]=min(mnmx_spont);

 p_v_width=(idx2-idx1)/32552*1000;  % DigiLynx sampling rate
else
 p_v_width=nan;
end
% Additional output
% if exist('lim1','var')
%     out.activation_start = lim1;
%     out.activation_end = lim2;
% end
end
% -------------------------------------------------------------------------
function mx = maxchannel(wv)

% Find largest channel
mean_wv = squeeze(nanmean(wv,1));   % mean waveform
amx = max(max(mean_wv));     % absolut maximum of mean waveforms
[mx,my] = find(mean_wv==amx,1,'first');      % mx: largest channel
end
