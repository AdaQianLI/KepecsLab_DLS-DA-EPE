function [H,F,J]= plot_mclust_projections2(cellid,selts_evoked,varargin)
%PLOT_MCLUST_PROJECTIONS2   Plot feature data for clusters.
%   PLOT_MCLUST_PROJECTIONS2(CELLID) plots feature data for a tagged cell
%   (CELLID) and it's tetrode pairs in all projections. Tagged cluster is
%   shown in orange and light-evoked spikes are overlayed in blue. Time
%   period for light-evoked activity is selected automatically (see
%   FINDSTIMPERIOD). Only spikes from the beginning of the first to the end
%   of the last stimulation protocol are included. These default behaviors
%   can be modified using the following optional input arguments
%   (parameter, value pairs, with default values):
%       'stim_period', [] - start and end of stimulation period after each
%           light pulse
%       'feature_names', 'Energy' - features for which feature data are
%           plotted; character or cell array
%       'marker', '+' - marker for the scatter plots
%       'marker_size', 2 - marker size for the scatter plots
%       'usefastplot', true - use fast plotting method (downsample points
%           to plot only one point per pixel; appears the same); faster,
%           but zoom is not implemented (for saving in image formats, e.g.
%           pdf or jpg); if false, full data is plotted (for viewing or
%           saving fig format)
%       'stimonly', true - only spikes from the beginning of the first to
%           the end of the last stimulation protocol are selected for
%           plotting; if false, all spikes are included
%       'plotlightspikes', true - if true, light-evoked spikes are 
%           superimposed
%
%   HS = PLOT_MCLUST_PROJECTIONS2(CELLID,...) returns the handles of the
%   figures in HS struct. The fields of HS are named according to the
%   features and channels: HS.(XFeatureXChannel_YFeatureYChannel), e.g.
%   HS.Amplitude1_Energy4.
%
%   Example:
%   plot_mclust_projections2(cellid,'feature_names',{'Amplitude','Energy'},...
%        'stim_period',[0.002 0.004]);
%
%   See also PLOTWAVEFORMS.

%   Sachin Ranade & Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   08-May-2012

%   Edit log: SPR 12/28/11; BH 5/8/12

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParameter(prs,'feature_names',{'feature_Energy'},@(s)iscell(s)|ischar(s))  % names of MClust features to use
addParameter(prs,'event','PulseOn_chr2',@ischar)   % reference event
addParameter(prs,'marker','.')  % marker for the plots
addParameter(prs,'marker_size',0.5,@isnumeric)  % marker size for the plots
addParameter(prs,'usefastplot',false,@(s)islogical(s)|ismember(s,[0 1]))  % fast plotting (no zoom)
addParameter(prs,'stimonly',true,@(s)islogical(s)|ismember(s,[0 1]))  % restrict to period between first and last light pulse
addParameter(prs,'plotlightspikes',true,@(s)islogical(s)|ismember(s,[0 1]))  % plot light-evoked spikes
addParameter(prs,'limwindow',[-0.02 0.2],...  % [-0.005 0.04]
    @(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
parse(prs,cellid,varargin{:})
g = prs.Results;
if ischar(g.feature_names)
    g.feature_names = {g.feature_names};
end


% cellidt=cellid;
% cellidt(cellidt=='.') = '_';
% datapath = getpref('cellbase','datapath');
% [subjectID, sessionID] = cellid2tags(cellid);
% resdir = [datapath '/' subjectID '/' sessionID 'Summary/'];
% if ~isdir(resdir)
%     mkdir(resdir)
% end
% % global DATAPATH
% % fs=filesep;
% % DATAPATH = uigetdir;
% % resdir = [DATAPATH 'Taggedprop'];

% Load stim events to get Pulse onset times.
ST = loadcb(cellid,'Stimevents');
pon = ST.(g.event)(~isnan(ST.(g.event)));

% Load spikes from Ntt file.
Nttfn = cellid2fnames(cellid,'Ntt');
all_spikes = LoadTT_NeuralynxNT(Nttfn);
TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
all_spikes = all_spikes;
val_spk_i = [find(all_spikes >= pon(1),1,'first') ...
    find(all_spikes <= pon(end),1,'last')]; % consider spikes only within the stimulation protocol to account for drift
nspk = length(all_spikes);
spk = loadcb(cellid,'Spikes');  % single unit all spikes
[junk,junk2,tagged_cell_inx] = intersect(spk,all_spikes);  %#ok<*ASGLU> % get indices for the cell
if ~isequal(junk,spk)  % check if all files have appropriate time stamps
    error('plot_mclust_projections:SpikeTimeMismatch','Mismatch between saved spike times and Ntt time stamps.')
end
if g.stimonly   % restrict to stimulation epoch
    tagged_cell_inx = tagged_cell_inx(tagged_cell_inx>val_spk_i(1)&tagged_cell_inx<val_spk_i(2));
end
tagged_cell_i = zeros(nspk,1);
tagged_cell_i(tagged_cell_inx) = 1;

% Cells on same tetrode including cellid.
[NumCell,tetpartners] = tetrodepairs(cellid);
tag_cell = strmatch(cellid,tetpartners);

% Spikes from each cell have an index. Spikes from noise have index 0.
cell_i = zeros(nspk,1);
for iCell = 1:NumCell
    spk = loadcb(tetpartners(iCell),'Spikes'); % load spike times.
    [junk,junk2,cell_inx] = intersect(spk,all_spikes); % get indices for the cell
    if g.stimonly   % restrict to stimulation epoch
        cell_inx = cell_inx(cell_inx>val_spk_i(1)&cell_inx<val_spk_i(2));
    end
    cell_i(cell_inx) = iCell;
end

% Load feature data for tetrode.
[r,s,t] = cellid2tags(cellid);
for k = 1:length(g.feature_names)
    prop = [g.feature_names{k} '.fd'];
    propfn = [getpref('cellbase','cell_pattern') num2str(t) '_' prop];
    sessionpath = cellid2fnames(cellid,'sess');
    propfn_path = [sessionpath filesep 'FD'];   % where the feature file can be found
    if ~isdir(propfn_path)
        propfn_path = sessionpath;
    end
    propfn_path = fullfile(propfn_path,propfn);
    wf_prop = load(propfn_path,'-mat');
    FeatureData(k,:,:) = wf_prop.FeatureData; %#ok<AGROW>
end

% Light-evoked spikes
if g.plotlightspikes
    
    % Latency of stimulated spikes
if isnan(g.limwindow(1)) || isnan(g.limwindow(2))
    g.limwindow(1) = 0.001;   % no activation detected
    g.limwindow(2) = 0.01;
end
    
    % Evoked spikes
%     tsegs_evoked = rel2abstimes(cellid,g.limwindow,'stim','PulseOn');   % convert period to epochs relative to pulses
%     selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find putative stimualated spikes
    [junk,junk2,evoked_cell_inx] = intersect(selts_evoked,all_spikes); % get indices for light-evoked spikes
end

% Plot
cmp = hsv(NumCell) / 4 + 0.75;
pcmb = allcomb(1:size(FeatureData,1),1:size(FeatureData,3));
cmb = flipud(combnk(pcmb(:,1)*10+pcmb(:,2),2));
NumComb = size(cmb,1);
for k = 1:NumComb
    fst = [floor(cmb(k,1)/10) mod(cmb(k,1),10)];
    scnd = [floor(cmb(k,2)/10) mod(cmb(k,2),10)];
    xdata = squeeze(FeatureData(fst(1),cell_i==0,fst(2)));
    ydata = squeeze(FeatureData(scnd(1),cell_i==0,scnd(2)));
    len=round(length(xdata)/8);
    xdata=xdata(1:len);
    ydata=ydata(1:len);
    
    % Open figure
%     namestr = [g.feature_names{fst(1)} num2str(fst(2)) '_'...
%         g.feature_names{scnd(1)} num2str(scnd(2))];
%     HS.(namestr) = figure('Position',[624 126 1092 852]);
  subplot(2,3,k);
hold on
    axis([min(xdata) max(xdata) min(ydata) max(ydata)])
    xlabel([g.feature_names{fst(1)} ': ' num2str(fst(2))])
    ylabel([g.feature_names{scnd(1)} ': ' num2str(scnd(2))])
    
    % Plot noise spikes in grey 
    if g.usefastplot
     H{k}= fastplot(xdata,ydata,[0.8 0.8 0.8],g.marker,g.marker_size/2);   % grey, all spikes from multi units, same rec. site 
    else
     H{k}= slowplot(xdata,ydata,[0.8 0.8 0.8],g.marker,g.marker_size/2);
    end
    
%     % Plot all clusters in color
%     for iC = 1:NumCell
%         xdatai = squeeze(FeatureData(fst(1),cell_i==iC,fst(2)));
%         ydatai = squeeze(FeatureData(scnd(1),cell_i==iC,scnd(2)));
%         if g.usefastplot
%             fastplot(xdatai,ydatai,cmp(iC,:),g.marker,g.marker_size);
%         else
%             slowplot(xdatai,ydatai,cmp(iC,:),g.marker,g.marker_size);
%         end
%     end
    
    % Plot tagged cluster
    xdatai = squeeze(FeatureData(fst(1),tagged_cell_i==1,fst(2)));
    ydatai = squeeze(FeatureData(scnd(1),tagged_cell_i==1,scnd(2)));
    if g.usefastplot
       F{k}= fastplot(xdatai(1:2:end),ydatai(1:2:end),[255 204 0]/255,g.marker,g.marker_size*2);  %orange, all spikes from single unit being selected
    else
       F{k}=slowplot(xdatai(1:2:end),ydatai(1:2:end),[255 204 0]/255,g.marker,g.marker_size*2);
    end
    
    % Plot light-evoked spikes
    if g.plotlightspikes
        xdatai = squeeze(FeatureData(fst(1),evoked_cell_inx,fst(2)));
        ydatai = squeeze(FeatureData(scnd(1),evoked_cell_inx,scnd(2)));
        if g.usefastplot
        J{k}=fastplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size*2);  %blue, light evoked spikes
        else
        J{k}=slowplot(xdatai,ydatai,[0 153 255]/255,'.',g.marker_size*2);
        uistack(J{k},'top')
        end
    end
%         fnm = [resdir cellidt '_mlustProjction.fig'];   % save
%         saveas(HS.(namestr),fnm);
end
% -------------------------------------------------------------------------
function C = allcomb(A,B)

% Convert to columns
A = A(:);
B = B(:);

% Combinations
as = arrayfun(@(k)horzcat(repmat(A(k),length(B),1),B),1:length(A),'UniformOutput',false);
C = cell2mat(as');

% -------------------------------------------------------------------------
function h = fastplot(x,y,clr,mrk,mrks)

% Reduce number of points
old_units = get(gca,'Units');
set(gca,'Units','pixels')
pos = get(gca,'Position');
xpixels = pos(3) + 1;
ypixels = pos(4) + 1;

xl = xlim;
mnx = xl(1);
mxx = xl(2);
yl = ylim;
mny = yl(1);
mxy = yl(2);
x2 = round((x-mnx)/(mxx-mnx)*(xpixels-1)) + 1;
y2 = round((y-mny)/(mxy-mny)*(ypixels-1)) + 1;
u = unique(x2*100000+y2);
y3 = mod(u,100000);
x3 = (u - y3) / 100000;
x4 = (x3 / xpixels) * (mxx - mnx) + mnx;
y4 = (y3 / ypixels) * (mxy - mny) + mny;

% Plot
h = plot(x4,y4,'.','Color',clr,'Marker',mrk,'MarkerSize',mrks);

% Restore axis units property
set(gca,'Unit',old_units)

% -------------------------------------------------------------------------
function h = slowplot(x,y,clr,mrk,mrks)

h = plot(x,y,'.','Color',clr,'Marker',mrk,'MarkerSize',mrks);