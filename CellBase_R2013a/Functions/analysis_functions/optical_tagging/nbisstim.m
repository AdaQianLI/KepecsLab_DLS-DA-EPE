function [p_value_act,Idiff_act,p_value_inh,Idiff_inh] = nbisstim(cellid)
%NBISSTIM   Assessment of optic tagging.
%   [P I] = NBISSTIM(CELLID) calculates information distances and
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
%   Endres and Schindelin, 2003)
%
%   NBISSTIM takes a different bin raster (typically aligned to 'BurstOn'
%   events) for baseline distribution than the the one for testing against
%   baseline (typically aligned to 'PulseOn' events). Number of test trials
%   used is maximized to 5000.
%
%   Reference:
%   Endres DM, Schindelin JE (2003) A new metric for probability
%   distributions. IEEE Transactions on Information Theory 49:1858-1860.
%
%   See also TAGGING, STIMES2BINRASTER and JSDIV.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com

% Input argument check
if nargin < 2
    win = [-0.02 0.06];  % [-0.6 0.6] time window for bin raster (for Hyun: 1.2)
    dt = 0.001;   % resolution of bin raster in s (for Hyun: 0.0005)
    dsply = 0;   % repress display
end

% Set parameters and load CellBase variables
EventName1 = 'BurstOn';
EventName2 = 'PulseOn';
ST = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'StimEvents');
epoch_pos1 = findcellstr(ST.events(:,1),EventName1);
epoch_pos2 = findcellstr(ST.events(:,1),EventName2);
if epoch_pos1 == 0 || epoch_pos2 == 0
    error('Epoch name not found');
end
stimes1 = ST.event_stimes{epoch_pos1};
stimes2 = ST.event_stimes{epoch_pos2};
time = win(1):dt:win(end);
valid_trials1 = find(~isnan(TE.(EventName1)));
% minfreq = min([TE.BurstNPulse]);
% maxpow = max([TE.PulsePower]);
% inx = ~isnan(TE.(EventName2)) & TE.BurstNPulse==minfreq & TE.PulsePower==maxpow;
% valid_trials2 = find(inx);
valid_trials2 = find(~isnan(TE.(EventName2)));
lm = 5000;    % downsaple if more pulses than 5000
if length(valid_trials2) > lm
    rp = randperm(length(valid_trials2));
    valid_trials2 = valid_trials2(sort(rp(1:lm)));
end

% Calculate bin rasters
spt1 = stimes2binraster(stimes1(valid_trials1),time,dt);
spt2 = stimes2binraster(stimes2(valid_trials2),time,dt);
sum_spt2=sum(spt2(:,21:end));
spt_baseline=spt2;
spt_mean=mean(sum_spt2);
spt_max=find(sum_spt2==max(sum_spt2));
spt_min=find(sum_spt2==min(sum_spt2));
if spt_max(1)<spt_min(1)
    spt_ref=spt_max(1);
else
    spt_ref=spt_min(1);
end
spt_baseline(:,21:spt_ref+5+21)=[];
%spt_baseline(:,end)=[];
[row_spt,col_spt]=size(spt_baseline);
% spt_baseline=reshape(spt_baseline,row_spt/5,col_spt*5);

% % % for i=1:200
% % %     B=spt_baseline(:,randperm(col_spt));
% % % spt_baseline_temp{1,i}=B;
% % % end
% % % [row_temp col_temp]=size(spt_baseline_temp);
% % % [row_ori col_ori]=size(spt_baseline);
% % % for i=1:200
% % %     spt_comp(:,(1+col_ori*(200-i)):col_ori*(201-i))=spt_baseline_temp{1,i};
% % % end

for i=1:100
    spt_baseline_temp{1,i}=spt_baseline(randperm(row_spt),randperm(col_spt));
end
[row_temp,col_temp]=size(spt_baseline_temp);
[row_ori,col_ori]=size(spt_baseline);
for i=1:100
    spt_comp(1:100,(1+col_ori*(100-i)):col_ori*(101-i))=spt_baseline_temp{1,i}(1:100,:);
end
spt_evoked=spt2(1:100,1:spt_ref+21);
% spt_evoked=spt2(1:100,21:spt_ref+5+20);
[row_evoked,col_evoked]=size(spt_evoked);
[row_comp,col_comp]=size(spt_comp);
col_comp2=floor(col_comp/col_evoked)*col_evoked;
spt_comp=spt_comp(:,1:col_comp2);
window=col_evoked*dt;
[p_value_act,Idiff_act] = salt2(spt_comp,spt_evoked,dt,window);
[p_value_inh,Idiff_inh] = salt2_inh(spt_comp,spt_evoked,dt,window);


% Set input arguments for rater plot and PSTH
if dsply
    SEvent1 = EventName1;
    SEvent2 = EventName2;
    FNum = 2;
    parts = 'all';
    sigma = 0.001;
    PSTHstd = 'on';
    ShEvent1 = {{EventName1}};
    ShEvent2 = {{EventName2}};
    ShEvColors = hsv(length(ShEvent1{1}));
    ShEvColors = mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
    
    % Plot raster plot and PSTH for 'BurstOn'
    figure
     set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName1,'SortEvent',SEvent1,'ShowEvents',ShEvent1,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on');
    pause(0.05);   % if reset the renderer two early, the same error occurs
   set(gcf,'renderer','opengl')   % reset renderer
    
    % Plot raster plot and PSTH for 'PulseOn'
    figure
    set(gcf,'renderer','painters')   % temporaray change renderer because OpenGL locks the plot which result an error in legend layout handling
    viewcell2b(cellid,'TriggerName',EventName2,'SortEvent',SEvent2,'ShowEvents',ShEvent2,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','on');
    pause(0.05);   % if reset the renderer two early, the same error occurs
    set(gcf,'renderer','opengl')   % reset renderer
end

