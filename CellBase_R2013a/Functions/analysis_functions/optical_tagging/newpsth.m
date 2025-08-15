function H = newpsth(cellid,spt,varargin)
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'win',[-0.1 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'event','PulseOn_chr2',@ischar)   % default reference event: 'PulseOn'\
addParameter(prs,'sigma',0.001,@isnumeric)
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds

parse(prs,cellid,varargin{:})
g = prs.Results;

margin = g.sigma * 3;     % add an extra margin to the windows


% Trial number and epoch length
[tno_baseline,tl] = size(spt);
[tno_test,tl] = size(spt);

% Pre-stimulus time window to consider for null hypothesis
stm = round(abs(g.win(1)-margin) / g.dt);

% Merged spike train
sptb = spt(:,1:stm);
[x0,allspks_baseline] = find(sptb);
ts_baseline = sort(allspks_baseline)';

sptt = spt(:,stm+1:end);
[x0,allspks_test] = find(sptt);
ts_test = stm + sort(allspks_test)';
ts = [ts_baseline ts_test];

% Calculate adaptive SDF with variable Gaussian Kernel
prob = [sum(sptb)/tno_baseline sum(sptt)/tno_test]  / g.dt;  % prob. if 1 ms bins; dt is measured in ms!
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
psth_aconv = agvd / g.dt * 1000;   % SDF
% psth_aconv = ppsth_aconv((stm+1+round((g.win(1)-margin)/g.dt)):(stm+1+fix((g.win(2)+margin)/g.dt)));
% psth_aconv_ave=normr(psth_aconv);
% prob = prob((stm+1+win(1)/g.dt):(stm+1+win(2)/g.dt));
% stn = abs(win(1)) / g.dt;

% Plot SDF
time = (g.win(1)-margin):g.dt:(g.win(2)+margin);
    subplot(2,1,1);
    cum_psth = cumsum(psth_aconv/sum(psth_aconv));
    H(1)=plot(time,cum_psth);
    xlim([time(1) time(end)]);
    ylim([0,1]);
    subplot(2,1,2);
    H(2)=cdfplot(psth_aconv);

end