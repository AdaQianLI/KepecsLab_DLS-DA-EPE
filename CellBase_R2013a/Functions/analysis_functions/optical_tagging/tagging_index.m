function  [p_value_act,Idiff_act,p_value_inh,Idiff_inh] = tagging_index(cellid,varargin)
%TAGGING_INDEX   Assessment of optic tagging.
%   [P I] = TAGGING_INDEX(CELLID) calculates information distances and
%   corresponding p values for light tagging for the cell given in CELLID
%   (see CellBase documentation). Briefly, a spike raster is calculated
%   with 1 ms resolution. The raster is devided to 10 ms time bins and
%   spike latency distribution for first spikes is computed within each
%   bin. Pairwise information divergence measures are calculated for the
%   before-light distributions to form a null-hypothesis distribution for
%   distances. The distances of the first after-light distribution (using
%   'BurstOn' events - see CellBase documentation) from all before-light
%   distributions is calculated and the median of these values is tested
%   against the null-hypothesis distribution. Output arguments P and I
%   correspond to a modified version of Jensen-Shannon divergence (see
%   Endres and Schindelin, 2003).
%
%   TAGGING_INDEX takes a different bin raster (typically aligned to
%   'BurstOn' events) for baseline distribution than the the one for
%   testing against baseline (typically aligned to 'PulseOn' events).
%   Number of test trials used is maximized to 5000 by default.
%
%   Default behavior of TAGGING_INDEX can be modified by using a set of
%   paramter, value pairs as optional input parameters. The folowing
%   parameters are implemented (with default values):
%       'event', 'PulseOn' - the event to which the window is locked
%   	'window', [-0.6 0.6] - extent of baseline and test period
%           relative to the event; in seconds
%   	'dt', 0.001 - time resolution of the bin raster; in seconds
%       'display', false - control of plotting event-locked raster plot
%       'event_filter', 'none' - filter light-stimulation trials; see
%           implemented filter types below
%       'maxtrialno', 5000 - maximal number of light-stimulation trials
%           included; if ther are more valid trials, they are randomly
%           down-sampled
%
%   Implemented event filters for 'stim' events:
%       'BurstNPulse' - only light pulse trains with a specified number of
%           pulses are used; this number has to be specified as a
%           'BurstNPulse', N (numeric) parameter, value pair given as a
%           struct (s.BurstNPulse = 1)
%       'BurstNPulse_maxPower' - the same as 'BurstNPulse' but light bursts
%           are further restricted to pulses with maximal power applied
%       'minNPulse_maxPower' - the light bursts with minimal number of
%           pulses and maximal power are used
%
%   Reference:
%   Endres DM, Schindelin JE (2003) A new metric for probability
%   distributions. IEEE Transactions on Information Theory 49:1858-1860.
%
%   See also TAGGING, STIMES2BINRASTER, FILTERTRIALS and JSDIV.

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'window_act',[-0.02 0.04],@(s)isnumeric(s)&isequal(length(s),2))  %[-0.6 0.6], time window for bin raster relative to the event, in seconds
addParameter(prs,'window_inh',[-0.02 0.01],@(s)isnumeric(s)&isequal(length(s),2))  %[-0.6 0.6], time window for bin raster relative to the event, in seconds
addParameter(prs,'binforsaltact',0.02,@isnumeric)  %[-0.6 0.6], time window for bin raster relative to the event, in seconds
addParameter(prs,'binforsaltinh',0.05,@isnumeric)  %[-0.6 0.6], time window for bin raster relative to the event, in seconds

addParameter(prs,'dt1',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParameter(prs,'dt2',0.001,@isnumeric)   % time resolution of the binraster, in seconds

addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
addParameter(prs,'event','PulseOn_chr2',@ischar)   % default reference event: 'PulseOn'\
addParameter(prs,'event_filter','none',@ischar)   % filter events based on properties
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
parse(prs,cellid,varargin{:})
g = prs.Results;

% g.event=TrigName1;
% g.window_act=actwin;
% g.window_inh=inhwin;
% g.binforsaltact=binforsaltact;
% g.binforsaltinh=binforsaltinh;
% g.dt1=dt1;
% g.dt2=dt2;

ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');

% % Spike times for baseline period
% if strcmpi(cellid,'default')
%     out = [strcat({'p_value_act','Idiff_act','p_value_inh','Idiff_inh'})];
%     return;
% end


try
epoch_pos1 = findcellstr(ST.events(:,1),g.event); % chr2
catch
epoch_pos1 =NaN;
end

if ~isnan (epoch_pos1)

stimes1 = ST.event_stimes{epoch_pos1};  
time1 = g.window_act(1):g.dt1:g.window_act(end);

valid_trials1 = filterTrials(cellid,'event_type','stim','event',g.event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);

% Downsaple if too many pulses
lm = g.maxtrialno;
if length(valid_trials1) > lm
    rp = randperm(length(valid_trials1));
    valid_trials1 = valid_trials1(sort(rp(1:lm)));
end

% Calculate bin rasters for
spt1 = stimes2binraster(stimes1(valid_trials1),time1,g.dt1); % ok<FNDSB>

% sum_spt1=sum(spt1);
% peak_idx=find(sum_spt1==max(sum_spt1));
baseline1_idx=(0-g.window_act(1))/g.dt1;
baseline1=spt1(:,1:baseline1_idx);
% spt_baseline(:,21:30)=[];
% spt_baseline(:,end)=[];
% spt_baseline=spt_baseline(1:row_idx,:);
[~,col_spt]=size(baseline1);
baseline1_=[];
% spt_baseline=spt_baseline(1:row_idx,:);

for K=1:1000
    a=randperm(col_spt,1);
    baseline1_=[baseline1_ baseline1(:,a)];
end
% [row_spt1,col_spt1]=size(baseline1);

evoked_idx1=(0-g.window_act(1))/g.dt1+1;
evoked_idx2=(g.window_act(2)-g.window_act(1))/g.dt1;
% evoked=spt1(1:row_spt1,evoked_idx1:evoked_idx2);
evoked1=spt1(:,evoked_idx1:evoked_idx2);


[p_value_act,Idiff_act] = salt2(baseline1_,evoked1,g.dt1,g.binforsaltact);


stimes2 = ST.event_stimes{epoch_pos1};  
time2 = g.window_inh(1):g.dt2:g.window_inh(end);

valid_trials2 = filterTrials(cellid,'event_type','stim','event',g.event,...
    'event_filter',g.event_filter,'filterinput',g.filterinput);

% Downsaple if too many pulses
lm = g.maxtrialno;
if length(valid_trials2) > lm
    rp = randperm(length(valid_trials2));
    valid_trials2 = valid_trials2(sort(rp(1:lm)));
end

% Calculate bin rasters for
spt2 = stimes2binraster(stimes2(valid_trials2),time2,g.dt2); % ok<FNDSB>

%%



% sum_spt1=sum(spt1);
% peak_idx=find(sum_spt1==max(sum_spt1));
baseline2_idx=(0-g.window_inh(1))/g.dt2;
baseline2=spt2(:,1:baseline2_idx);
% spt_baseline(:,21:30)=[];
% spt_baseline(:,end)=[];
[~,col_spt]=size(baseline2);
baseline2_=[];
% spt_baseline=spt_baseline(1:row_idx,:);
for K=1:1000
    a=randperm(col_spt,1);   
    baseline2_=[baseline2_ baseline2(:,a)];
end

evoked_idx1=(0-g.window_inh(1))/g.dt2+1;
evoked_idx2=(g.window_inh(2)-g.window_inh(1))/g.dt2;
evoked2=spt2(:,evoked_idx1:evoked_idx2);

[p_value_inh,Idiff_inh] = salt2_inh(baseline2_,evoked2,g.dt2,g.binforsaltinh);

% Set input arguments for rater plot and PSTH
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
%     % Plot raster plot and PSTH for 'BurstOn'
% % %     figure
% % %     set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
% % %     viewcell2b(cellid,'TriggerName',EventName1,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
% % %         'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
% % %         'EventMarkerWidth',0,'PlotZeroLine','on')
% % %     pause(0.05)   % if reset the renderer two early, the same error occurs
% % %     set(gcf,'renderer','opengl')   % reset renderer
% % %     
% % %     % Plot raster plot and PSTH for 'PulseOn'
% % %     figure
% % %     set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
% % %     viewcell2b(cellid,'TriggerName',EventName2,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
% % %         'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
% % %         'EventMarkerWidth',0,'PlotZeroLine','on')
% % %     pause(0.05)   % if reset the renderer two early, the same error occurs
% % %     set(gcf,'renderer','opengl')   % reset renderer
% end

else
    p_value_act=NaN;
    Idiff_act=NaN;    
    p_value_inh=NaN;
    Idiff_inh=NaN;
end
end

% out=[p_value_act,Idiff_act,p_value_inh,Idiff_inh];
% % % Calculate information distances and p values
% % res = 10;   % resolution in ms
% % dtt = g.dt * 1000;   % resolution of bin raster in ms
% % wn = g.window * 1000;   % window boundaries in ms
% % [p_value Idiff] = isikldist(spt1,spt2,dtt,wn,res);

% -------------------------------------------------------------------------
% function [p_value Idiff] = isikldist(spt_baseline,spt_test,dt,win,res)
% 
% % Trial number and epoch length
% [tno tl] = size(spt_baseline); %#ok<NASGU>
% 
% % Number of bins for ISI histograms
% nmbn = round(res/dt);
% 
% % Pre-stimulus time window to consider for null hypothesis
% st = abs(win(1)) / dt;   % number of pre-stim values in 'spt'
% 
% % ISI histogram - baseline
% edges = 0:nmbn+1;
% nm = floor(st/nmbn);
% lsi = zeros(tno,nm);   % ISI's
% slsi = zeros(tno,nm);  % sorted ISI's
% hlsi = zeros(nmbn+1,nm);    % ISI hist.; +1: zero when no spike in the segment
% nhlsi = zeros(nmbn+1,nm);   % normalized ISI histogram 
% next = 1;
% for t = 1:nmbn:st
%     for k = 1:tno
%         cspt = spt_baseline(k,t:t+nmbn-1);
%         pki = find(cspt,1,'first');
%         if ~isempty(pki)
%             lsi(k,next) = pki;
%         else
%             lsi(k,next) = 0;
%         end
%     end
%     slsi(:,next) = sort(lsi(:,next));
%     hst = hist(slsi(:,next),edges);
%     hlsi(:,next) = hst(1:end-1);
%     nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));
%     next = next + 1;
% end
% 
% % ISI histogram - test
% tno_test = size(spt_test,1);
% lsi_tt = nan(tno_test,1);
% for k = 1:tno_test
%     cspt = spt_test(k,st+1:st+nmbn);
%     pki = find(cspt,1,'first');
%     if ~isempty(pki)
%         lsi_tt(k,1) = pki;
%     else
%         lsi_tt(k,1) = 0;
%     end
% end
% slsi_tt = sort(lsi_tt(:,1));
% hst = hist(slsi_tt,edges);
% hlsi(:,next) = hst(1:end-1);
% nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));
% 
% % figure      % plot ISIs
% % imagesc(lsi)
% % figure      % plot sorted ISIs
% % imagesc(slsi)
% % figure      % plot ISI histograms
% % imagesc(hlsi(2:end,:))
% 
% % Symmetric KL-divergence and JS-divergence
% kn = st / nmbn + 1;
% jsd = nan(kn,kn);  % pairwise modified JS-divergence (which is a metric!)
% for k1 = 1:kn
%     D1 = nhlsi(:,k1);
%     for k2 = k1+1:kn
%         D2 = nhlsi(:,k2);
%         jsd(k1,k2) = sqrt(JSdiv(D1,D2)*2);
%     end
% end
% % figure    % plot KL-distance
% % imagesc(kld)
% 
% % Calculate p-value and information difference
% [p_value Idiff] = makep(jsd,kn);
% % keyboard
% 
% % -------------------------------------------------------------------------
% function [p_value Idiff] = makep(kld,kn)
% % Calculates p value from distance matrix.
% 
% pnhk = kld(1:kn-1,1:kn-1);
% nullhypkld = pnhk(~isnan(pnhk));   % nullhypothesis
% testkld = median(kld(1:kn-1,kn));  % value to test
% sno = length(nullhypkld(:));   % sample size for nullhyp. distribution
% p_value = length(find(nullhypkld>=testkld)) / sno;
% Idiff = testkld - median(nullhypkld);