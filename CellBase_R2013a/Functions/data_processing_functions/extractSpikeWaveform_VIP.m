function wave= extractSpikeWaveform_VIP(cellid,SpikeTimes,varargin)
%EXTRACTSPIKEWAVEFORMS   Peak-aligned waveforms.
%   WAVE = EXTRACTSPIKEWAVEFORMS(CELLID,SPIKETIMES) returns waveforms for 
%   selected spikes (SPIKETIMES) of a given cell (CELLID). Spikes are
%   algined based on peak times of the tetrode channel with largest mean
%   amplitude.
%
%   WAVE = EXTRACTSPIKEWAVEFORMS(CELLID,SPIKETIMES,PARAM,VALUE) accepts
%   additional optional input arguments in the form of parameter, value
%   pairs (with default values):
%       'chans', 'all' - filter for tetrode channels (see below details of
%           implmented filter types
%       'smooth', 'none' - controls smoothing of the waveforms; implemented
%           only for average waveforms; default: no smoothing; alternative:
%           'spline': spline interpolation
%       'usemedian', false - if true, median is used for averaging
%           waveforms instead of mean
%
%   Filter types (possible values for 'chans'):
%       'all' - all waveforms
%       'max' - all waveforms on largest channel
%       'min' - all waveforms on smallest channel
%       'mean_all' - mean/median waveform on all channels
%       'mean_max' - mean/median waveform on largest channel
%
%   See also PLOTWAVEFORMS.

%   Edit log: BH 4/26/12

% Input arguments
default_args={...
    'chans',           'all';...     % filter tetrode channels
    'smooth',          'none';...  % control whether to smooth; implemented only for average outputs
    'usemedian',       false;... % if output is smoothed, either mean or median is applied
    };
[g, error] = parse_args(default_args,varargin{:});

% Load waveform data (Ntt file)
Datfile = cellid2fnames(cellid,'dat');
% TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
[all_spikes,all_waves] = mClustTrodes4LoadingEnginey(Datfile);
%%[junk junk2 evoked_inx] = intersect(SpikeTimes,all_spikes*TIMEFACTOR);
[junk,junk2,evoked_inx] = intersect(SpikeTimes,all_spikes);
if ~isequal(junk,SpikeTimes)   % internal check for spike times
 %%error('extractSpikeWaveforms:SpikeTimeMismatch','Mismatch between extracted spike times and Ntt time stamps.')
SpikeTimes_matched=SpikeTimes(junk2);
else
  SpikeTimes_matched=SpikeTimes;  
end

sel_wv = all_waves(evoked_inx,:,:);   % selected waveforms

% Align waveforms
aligned_wv = wvalign(sel_wv);

% Apply filter on the tetrode channels
if g.usemedian
    avegageFcn = @nanmedian;    % define which function to use for averaging
else
    averageFcn = @nanmean;
end
if isnumeric(g.chans),
    wave = squeeze(aligned_wv(:,g.chans,:));    % return specified channels
    g.smooth = 'none';   % smoothing is only implemented for average outputs
else
    switch g.chans
        case 'all'
            wave = aligned_wv;   % all waveforms
            g.smooth = 'none';   % smoothing is only implemented for average outputs
        case 'max'
            mx = maxchannel(aligned_wv);     % mx: largest channel
            wave = squeeze(aligned_wv(:,mx,:));     % all waveforms on largest channel
            g.smooth = 'none';   % smoothing is only implemented for average outputs
        case 'min'
            mn = minchannel(aligned_wv);     % mn: smallest channel
            wave = squeeze(aligned_wv(:,mn,:));     % all waveforms on smallest channel
            g.smooth = 'none';   % smoothing is only implemented for average outputs
        case 'mean_all'
            wave = squeeze(averageFcn(aligned_wv,1));   % mean/median waveform on all channels
        case 'mean_max'
            mx = maxchannel(aligned_wv);     % mx: largest channel
            wave = averageFcn(squeeze(aligned_wv(:,mx,:)));    % mean/median waveform on largest channel
    end
end

% Smooth waveforms
if isequal(g.smooth,'spline')
    wave = wvsmooth(wave);
end
end
% -------------------------------------------------------------------------
function aligned_wv = wvalign(wv)

% Find largest channel
mx = maxchannel(wv);     % mx: largest channel

% Align spikes to peak
n = size(wv,1);    % number of spikes
aligned_wv = nan(n,4,96);   % extended waveforms
for k = 1:n
     [mval,minx] = max(wv(k,mx,:));
     shift = 40 - minx;   % shift is decided based on the largest channel
     aligned_wv(k,mx,shift:shift+31) = wv(k,mx,:);
     for oc = setdiff(1:4,mx)   % align other channels with the same shift
         aligned_wv(k,oc,shift:shift+31) = wv(k,oc,:);
     end
end
aligned_wv = aligned_wv(:,:,32:63);     % transform back to original size
end
% -------------------------------------------------------------------------
function smoothed_wv = wvsmooth(wv)

% Smooth with spline
x = 0:31;
xx = 0:.1:31;
smoothed_wv = spline(x,wv,xx);
end
% -------------------------------------------------------------------------
function mx = maxchannel(wv)

% Find largest channel
mean_wv = squeeze(nanmean(wv,1));   % mean waveform
amx = max(max(mean_wv));     % absolut maximum of mean waveforms
[mx my] = find(mean_wv==amx,1,'first');     % mx: largest channel
end

% -------------------------------------------------------------------------
function mn = minchannel(wv)

% Find smallest channel
mean_wv = squeeze(nanmean(wv,1));   % mean waveform
[mnv mn] = min(max(mean_wv,[],2));     % mnv: minimal peak value; mn: smallest channel    
end