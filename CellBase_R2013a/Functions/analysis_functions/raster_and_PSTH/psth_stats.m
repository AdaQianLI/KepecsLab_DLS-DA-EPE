function [stats,s,y]= psth_stats(cellid,spt,psth,varargin)
%PSTH_STATS   Statistics for firing rate change.
%   PSTH_STATS tests whether spike number after an event changes.
%
%   STATS = PSTH_STATS(spt,PSTH,DT,WIN) finds minimal/maximal firing as
%   minimum/maximum PSTH within 100 ms from time zero (determined by time
%   limits of the window around 0, WIN and temporal resolution, DT).
%   Baseline firing is determined by mean firing probability from -250 ms 
%   to 0 using spt, the bin raster corrsponding to the PSTH. Next, the time
%   course of inhibition/activation is assessed by crossings of the
%   half-distence between the extreme and the baseline before and after the
%   minimum/maximum. This temporal window of inhibition/activation is then
%   used to find corresponding intervals around local extremes in the
%   baseline raster. Spike counts for baseline and spike counts in the
%   previously determined inhibition window are compared using Mann-Whitney
%   U-test (one-sided p-value is calculated; multiply it by 2 for two-sided
%   p-value).
%
%   Default behavior of PSTH_STATS can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'baselinewin', [-0.25 0] - limits of baseline window for 
%           statistical testing, time relative to 0 in seconds
%   	'testwin', [0 0.1] - limits of test window for statistical testing,
%           time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals; in proportion of
%           the peak-baseline difference
%       'display', false - controls plotting.
%
%   The output structure STATS contains the following fields:
%       'baseline' - baseline firing probability
%       'minvalue' - minimal firing rate in the test window
%       'inhibition_start' - start time of inhibition in seconds
%       'inhibition_end' - end time of inhibition in seconds
%       'inhibition_peak' - peak time of inhibition in seconds
%       'inhibition_time' - duration of inhibition in seconds
%       'Wpi' - p value for Mann-Whitney test for significant inhibition
%       'maxvalue' - maximal firing rate in the test window
%       'activation_start' - start time of activation in seconds
%       'activation_end' - end time of activation in seconds
%       'activation_peak' - peak time of activation in seconds
%       'activation_time' - duration of activation in seconds
%       'Wpa' - p value for Mann-Whitney test for significant activation
%
%   IMPORTANT NOTE: Please note that the baseline window is split to
%   smaller windows of the size of the test window. Use baseline windows
%   that are somewhat bigger than an integer multiple of the test window
%   length. The baseline window will be cropped BEFORE the nearest integer
%   multiple of the test window size. For instance, a test window of [0
%   0.1] and a baseline window of [-0.22 0] will result in an effective
%   baseline window of [-0.2 0], corresponding to twice the size of the
%   baseline window.
%
%   See also ULTIMATE_PSTH.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   12-Aug-2012

%   Edit log: BH 8/12/12, 8/27/12, 7/15/13

% Default arguments
prs = inputParser;
% addRequired(prs,'spt')   % bin raster
% addRequired(prs,'psth')   % PSTH
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster and PSTH, in seconds
addParameter(prs,'sigma',0.001,@isnumeric)
addParameter(prs,'win',[-0.1 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'plotwin',[-0.1 0.1],@(s)isnumeric(s)&isequal(length(s),2))
addParameter(prs,'baselinewin',[-0.25 0],@(s)isnumeric(s)&isequal(length(s),2))  % [-0.25 0],baseline time window relative to the event for stat. testing, in seconds
addParameter(prs,'testwin',[0 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % test time window relative to the event for stat. testing, in seconds
addParameter(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>0&s<1)   % 0.5, threshold used to assess interval limits
addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
addParameter(prs,'FigureNum',1,@isnumeric) 



parse(prs,cellid,varargin{:})
g = prs.Results;

margin = g.sigma * 3;     % add an extra margin to the windows

% g.testwin(2)=g.testwin(2)+margin;
g.baselinewin(1)=g.baselinewin(1)-margin;
g.win(1)=g.win(1)-margin;
g.win(2)=g.win(2)+margin;

xscale=g.plotwin;
g.plotwin(1)=g.plotwin(1)-margin;
g.plotwin(2)=g.plotwin(2)+margin;

% Error handling
twlen = diff(g.testwin);   % length of the test window
bwlen = diff(g.baselinewin);   % length of the baseline window
if bwlen < twlen
    warning('psth_stats:inputArg','PSTH_STATS: For a valid MW-test, the baseline window should be longer than the test window.')
end
if bwlen > -g.win(1)
    error('psth_stats:inputArg','PSTH_STATS: The PSTH window does not cover the baseline window.')
end
if twlen > g.win(2)
    error('psth_stats:inputArg','PSTH_STATS: The PSTH window does not cover the test window.')
end

% Index for time 0
st = fix(abs(g.win(1)) / g.dt);   % in ms
nullindex = st + 1;
startidx=fix((g.plotwin(1)-g.win(1))/g.dt)+1;
endidx=fix(g.plotwin(2)/g.dt)+nullindex;

% Window for testing the potential effect
WNb = [fix(g.baselinewin(1)/g.dt)+nullindex fix(g.baselinewin(2)/g.dt)+nullindex-1];   % baseline window; convert to indices; exlude 0
WNt = [fix(g.testwin(1)/g.dt)+nullindex fix(g.testwin(2)/g.dt)+nullindex];   % test window; convert to indices
% WNb = fix(WNb)+1;
% WNt = fix(WNt);
lWNb = WNb(2) - WNb(1) + 1;   % length of baseline window
lWNt = WNt(2) - WNt(1) + 1;   % length of test window

% Trial number
tno = size(spt,1);

% Time vector
time = g.win(1):g.dt:g.win(2);

% Spiking probability
sptb = spt(:,1:st);    % st < time 0, baseline spikes
sptt = spt(:,st+1:end);   % st+1 >time 0, test spikes
% probb = sum(sptb) / tno;   % spiking prob. (strictly) before time 0

% Inhibition time
% psth_baseline = mean(probb(WNb(1):WNb(2))) / g.dt;  % spikes/sec (was spikes/bin before)
psth_baseline = mean(psth(WNb(1):WNb(2)));  
psth_minafter = min(psth(WNt(1):WNt(2)));
% if psth_minafter > psth_baseline     % putative inhibition, if it goes below baseline
%     inhibition_start = NaN;
%     inhibition_end = NaN;
%     inhibition_peak = NaN;
%     inhibition_time = 0;   % if firing does not go below baseline
%     Wpi = NaN;
%     Whi=NaN
% else
    mininx = WNt(1)-1+find(psth(WNt(1):WNt(2))==psth_minafter,1,'first');   % minimal firing
    thr = psth_baseline - (psth_baseline - psth_minafter) * g.relative_threshold;  % threshold is determined in proportion of peak-baseline distance
    pis = valuecrossing(time(WNt(1):mininx),psth(WNt(1):mininx),thr,'down');
    pis_inx = valuecrossing(WNt(1):mininx,psth(WNt(1):mininx),thr,'down');
    if isempty(pis)
        pis = time(WNt(1));
        pis_inx = WNt(1);
    end
    pis_inx = fix(pis_inx(end));
    inhibition_start = pis(end);   % last crossing of half-baseline probability before minimum
    pie = valuecrossing(time(mininx:WNt(2)),psth(mininx:WNt(2)),thr,'up');
    pie_inx = valuecrossing(mininx:WNt(2),psth(mininx:WNt(2)),thr,'up');
    if isempty(pie)
        pie = time(WNt(2));
        pie_inx = WNt(2);
    end
    pie_inx = fix(pie_inx(1));
    inhibition_end = pie(1);   % first crossing of half-baseline probability after minimum
    inhibition_time = inhibition_end - inhibition_start;
    inhibition_peak = time(mininx) - time(st+1);    % peak time of inhibition
    
    % Nullhypothesis distribution
    wns = pie_inx - pis_inx + 1;
    wnnm = floor(lWNb/lWNt);   % split up the baseline according to test window length
    psp = nan(wns,wnnm*tno);
    for k = 1:wnnm
        inx = WNb(2)-k*lWNt+1:WNb(2)-(k-1)*lWNt;
        cwn = sptb(:,inx);
        cwnps = psth(inx);
        mcw = find(cwnps==min(cwnps));
        mcw = mcw(1);
        inx2 = mcw-floor(wns/2):mcw+ceil(wns/2)-1;
        if ~isempty(inx2)
        inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));       
        if ~ismember(mcw,inx2)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        mcwn = cwn(:,inx2);
        if ~isequal(size(mcwn,2),wns)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        psp(:,(k-1)*tno+1:k*tno) = mcwn';
        else
           psp=NaN; 
        end
    end
%     if any(isnan(psp))
%         error('vipisinfluenced2:nullhypoIndexing','Programming error.')
%     end
    spno_null = sum(psp);
    
    % Test distribution
    spno_test = sum(sptt(:,pis_inx-st:pie_inx-st),2);
    
    % Mann-Whitney test
    [Wpi,Whi] = b_ranksum2(spno_null,spno_test,'alpha',0.01);
    ranks = tiedrank([spno_null spno_test']);
    tp = length(spno_null);
    nullranksum = mean(ranks(1:tp));
    testranksum = mean(ranks(tp+1:end));
    if testranksum > nullranksum    % one-sided test
        Wpi = NaN;
        Whi = 0;
    end
    if Whi>0
        clri = [0 153 255] / 256;
    else
        clri = [102 255 255] / 256;
    end
% end

    
% Activation time
maxafter = max(psth(WNt(1):WNt(2)));
% if maxafter < psth_baseline     % putative activation, if firing goes above baseline
%     activation_start = NaN;
%     activation_end = NaN;
%     activation_peak = NaN;
%     activation_time = 0;   % if firing does not go above baseline
%     Wpa = NaN;
%     Wha = NaN;
% else
    maxinx = WNt(1) - 1+find(psth(WNt(1):WNt(2))==maxafter,1,'first');   % maximal firing
    thr = psth_baseline + (maxafter - psth_baseline) * g.relative_threshold;  % threshold is determined in proportion of peak-baseline distance
    pas = valuecrossing(time(WNt(1):maxinx),psth(WNt(1):maxinx),thr,'up');
    pas_inx = valuecrossing(WNt(1):maxinx,psth(WNt(1):maxinx),thr,'up');
    if isempty(pas)
        pas = time(WNt(1));
        pas_inx = WNt(1);
    end
    pas_inx = round(pas_inx(end));
    activation_start = pas(end);   % last crossing of one and a half-baseline probability before maximum
    pae = valuecrossing(time(maxinx:WNt(2)),psth(maxinx:WNt(2)),thr,'down');
    pae_inx = valuecrossing(maxinx:WNt(2),psth(maxinx:WNt(2)),thr,'down');
    if isempty(pae)
        pae = time(WNt(2));
        pae_inx = WNt(2);
    end
    pae_inx = round(pae_inx(1));
    activation_end = pae(1);   % first crossing of one and a half-baseline probability after maximum
    activation_time = activation_end - activation_start;
    activation_peak = time(maxinx) - time(st+1);    % peak time of activation
    
    % Nullhypothesis distribution
    wns = pae_inx - pas_inx + 1;
    wnnm = floor(lWNb/lWNt);   % split up the baseline according to test window length
    psp = nan(wns,wnnm*tno);
    for k = 1:wnnm
        inx = WNb(2)-k*lWNt+1:WNb(2)-(k-1)*lWNt;
        cwn = sptb(:,inx);
        cwnps = psth(inx);
        mcw = find(cwnps==max(cwnps));
        mcw = mcw(1);
        inx2 = mcw-floor(wns/2):mcw+ceil(wns/2)-1;
       if ~isempty(inx2)
           inx2 = inx2 - min(0,inx2(1)-1) - (max(length(inx),inx2(end)) - length(inx));
        if ~ismember(mcw,inx2)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        mcwn = cwn(:,inx2);
        if ~isequal(size(mcwn,2),wns)
            error('vipisinfluenced2:nullhypoIndexing','Programming error.')
        end
        psp(:,(k-1)*tno+1:k*tno) = mcwn';
       else
           psp=NaN; 
        end
    end
    if any(isnan(psp))
        error('vipisinfluenced2:nullhypoIndexing','Programming error.')
    end
    spno_null = sum(psp);
    
    % Test distribution
    spno_test = sum(sptt(:,pas_inx-st:pae_inx-st),2);
    
    % Mann-Whitney test
    [Wpa,Wha] = b_ranksum2(spno_null,spno_test,'alpha',0.01);
    ranks = tiedrank([spno_null spno_test']);
    tp = length(spno_null);
    nullranksum = mean(ranks(1:tp));
    testranksum = mean(ranks(tp+1:end));
    if testranksum < nullranksum    % one-sided test
        Wpa = NaN;
        Wha = 0;
    end
    if Wha>0
        clra = 'red';
    else
        clra = [255 102 0] / 256;
    end
% end

% Plot
if g.display

    plot(time(startidx:endidx),psth(startidx:endidx),'k');
    xlim(xscale);
    
    if exist('clri','var')
        hold on
        s(1)=plot(time(pis_inx:pie_inx),psth(pis_inx:pie_inx),'Color',clri,'LineWidth',2);
        x_lim = xlim;
        y_lim = ylim;
        y(1)=text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,['{\itMW test}, p = ',num2str(Wpi)],'Color',clri,'FontSize',8);
    else
        s(1)=plot(0,0,'.w');
         x_lim = xlim;
        y_lim = ylim;
        y(1)=text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,'NaN','FontSize',8);
    end
    if exist('clra','var')
        hold on
        s(2)=plot(time(pas_inx:pae_inx),psth(pas_inx:pae_inx),'Color',clra,'LineWidth',2);
        x_lim = xlim;
        y_lim = ylim;
        y(2)=text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,['{\itMW test}, p = ',num2str(Wpa)],'Color',clra,'FontSize',8);
    else
        s(2)=plot(0,0,'.w');
         x_lim = xlim;
        y_lim = ylim;
        y(2)=text(x_lim(1)+(x_lim(2)-x_lim(1))*0.6,y_lim(1)+(y_lim(2)-y_lim(1))*0.8,'NaN','FontSize',8);

    end
end

% Output statistics
stats.baseline = psth_baseline;
stats.minvalue = psth_minafter;
stats.inhibition_start = inhibition_start;
stats.inhibition_end = inhibition_end;
stats.inhibition_peak = inhibition_peak;
stats.inhibition_time = inhibition_time;
stats.Wpi = Wpi;
stats.maxvalue = maxafter;
stats.activation_start = activation_start;
stats.activation_end = activation_end;
stats.activation_peak = activation_peak;
stats.activation_time = activation_time;
stats.Wpa = Wpa;