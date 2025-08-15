function error_list = taggedprop(I,issave)
%TAGGEDPROP   Properties of putative tagged neurons.
%   TAGGEDPROP(I,ISSAVE) performs analysis aiding the definite decision
%   about tagging. Input parameters: 
%       I - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, putative tagged
%           cells are selected (ID>20, L-ratio<0.15, H index<0.01, R>0.9;
%           see TAGGING, LRATIO, NBISSTIM and SPIKESHAPECORR for
%           details on these measures)
%       ISSAVE - controls saving
%
%   The following analyses are performed (output variables saved in mat
%   files and figures in pdf):
%       Reliability, latency and jitter of evoked spikes after 'BurstOn',
%           'PulseOn' and frequency-restricted 'PulseOn' events; see
%           RELIABILITY_LATENCY_JITTER for details
%       H index for 'BurstOn', 'PulseOn' and frequency-restricted 'PulseOn'
%           events; see TAGGING_INDEX for details
%       Waveforms of spontaneous and light-evoked spikes; see PLOTWAVEFORMS
%           for details
%       All projections of feature data in the Energy-Amplitue space with
%           the putative tagged cells shown in orange and the light-evoked 
%           spikes overlayed in blue; see PLOT_MCLUST_PROJECTIONS2 for
%           details
%       Cluster quality indices restricted to light-evoked spikes in the
%           Energy-Amplitude as well as in the Energy-WavePC1 space; see
%           LIGHTCLUSTER and LRATIO for details
%       Raster plot and PSTH aligned to 'BurstOn' and 'PulseOn' events
%           (only 1000 pulses shown in the latter); see PLOT_RASTER_PSTH,
%           VIEWCELL2B and PLOT_RASTER2A for details.
%
%   ERROR_LIST = TAGGEDPROP(I,ISSAVE) returns a structure with all caught
%   errors. ERROR_LIST includes cell IDs, error messages and the captured
%   exception variables. A time stamped ERROR_LIST is also assigned in base
%   workspace and saved to the results directory automatically.
%
%   See also TAGGING, LRATIO, SPIKESHAPECORR, RELIABILITY_LATENCY_JITTER,
%   TAGGING_INDEX, PLOTWAVEFORMS, PLOT_MCLUST_PROJECTIONS2, LIGHTCLUSTER,
%   PLOT_RASTER_PSTH, VIEWCELL2B and PLOT_RASTER2A.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   05-Oct-2012

%   Edit log: BH 5/10/12, 4/19/13

% Pass the control to the user in case of error
dbstop if error
 
% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1
    I = [];
end
 
% Directories
global DATAPATH
DATAPATH = uigetdir;
resdir = [DATAPATH 'Summary'];
if ~isdir(resdir)
    mkdir(resdir)
end

% Load CellBase
load(getpref('cellbase','fname'),'CELLIDLIST');

    if isnumeric(I)
        putative_tagged = CELLIDLIST(I);   % index set to CELLIDLIST
    elseif ischar(I)
        putative_tagged = {I};   % only one cellID passed
    elseif iscellstr(I)
        putative_tagged = I;   % list of cell IDs
    else
        error('taggedprop:inputArg','Unsupported format for cell IDs.')
    end

% All the analysis promised in the help
NumCells = length(putative_tagged);
error_list = struct('cellid',{},'message',{},'exception',{});  % keep a list of caught errors 
errinx = 0;
for k = 1:NumCells
    cellid = putative_tagged{k};
    try
        main(cellid,issave,resdir)  % every analysis
    catch ME
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
        errinx = errinx + 1;  % error counter
        error_list(errinx).cellid = cellid;   % record problematic cell ID
        error_list(errinx).message = ME.message;   % error message
        error_list(errinx).exception = ME;   % exception structure
    end
end

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
list_name = ['error_list_' dsr];  % variable name
list_name = regexprep(list_name,'-','_');
list_name = regexprep(list_name,' ','_');
assignin('base',list_name,error_list)   % assign error list in base workspace
error_fname = fullfile(resdir,[list_name '.mat']);   % file name
save(error_fname,'error_list')   % save error list

% -------------------------------------------------------------------------
function main(cellid,issave,resdir)

% Add 'PulseOn' event if missing
% % ST = loadcb(cellid,'STIMSPIKES');
% % if isequal(findcellstr(ST.events(:,1),'PulseOn'),0)
% %     [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
% %     prealignSpikes_hp(cellid,'events',...
% %         {'PulseOn' 'PulseOn' 'PulseOn' [lim1 lim2]},...
% %         'epochs',[],'filetype','stim','ifsave',1,'ifappend',1)
% % else
% %     prealignSpikes_hp(cellid,'events',...
% %         {'PulseOn' 'PulseOn' 'PulseOn' [lim1 lim2]},'epochs',[],...
% %         'filetype','stim','ifsave',1,'writing_behavior','replace');
% % end

SE = loadcb(cellid,'StimEvents');
burst_types = sort(unique(SE.BurstNPulse),'ascend');
burst_types(isnan(burst_types)) = [];
NumBurstTypes = length(burst_types);
% H-index for 'PulseOn'
[Hindex_act D_KL_act Hindex_inh D_KL_inh] = tagging_index(cellid,'event','PulseOn');
% H-index for 'BurstOn'
[Hindex_act_burston D_KL_act_burston Hindex_inh_burston D_KL_inh_burston] = tagging_index(cellid,'event','BurstOn');  %#ok<NASGU> % H-index, D_KL
% Spike shape correlation
R = spikeshapecorr(cellid);
% Efficiency, latency and jitter for 'PulseOn'
[E_pulseon L_pulseon J_pulseon B_pulseon M_pulseon lim_act1 lim_act2 Wp_act MinRate lim_inh1 lim_inh2 Wp_inh]=reliability_latency_jitter(cellid,'event','PulseOn');
% efficiency,latency,jitter,baseline firing rate,max rate,           p for mann whitney test                                                                        
if (M_pulseon-B_pulseon)>(B_pulseon-MinRate)
    A1=lim_act1;
    A2=lim_act2;
else
    A1=lim_inh1;
    A2=lim_inh2;
end

% Efficiency, latency and jitter for 'BurstOn'
[E_burston L_burston J_burston B_burston M_burston] = reliability_latency_jitter(cellid,'event','BurstOn');

% Calculate the same tagging variables for bursts of different frequencies
Hindex_frequency = nan(1,NumBurstTypes);
D_KL_frequency = nan(1,NumBurstTypes);
E_frequency = nan(1,NumBurstTypes);
L_frequency = nan(1,NumBurstTypes);
J_frequency = nan(1,NumBurstTypes);
B_frequency = nan(1,NumBurstTypes);
M_frequency = nan(1,NumBurstTypes);
for bt = 1:NumBurstTypes
    
    fi = struct('BurstNPulse',burst_types(bt));
    [Hindex_act_frequency(bt) D_KL_act_frequency(bt) Hindex_inh_frequency(bt) D_KL_inh_frequency(bt)] = ...
        tagging_index(cellid,'event_filter','BurstNPulse_maxPower','filterinput',fi);  % H-index, D_KL
    [E_frequency(bt) L_frequency(bt) J_frequency(bt) B_frequency(bt) M_frequency(bt)] = reliability_latency_jitter(cellid,...
        'event_filter','BurstNPulse_maxPower','filterinput',fi);  % efficiency, latency, jitter
    
end

%%%spike waveform width
[~, waveprop,wave_spont wave_evoked]= plotwaveforms(cellid,'maxnum',5000); % 'maxnum',  2000 - maximum number of spikes to plot

save([resdir '\TAGGEDPROP_' regexprep(cellid,'\.','_') 'spkprop.mat'],...
        'waveprop','wave_spont','wave_evoked');
    
% % Plot efficiency
% HE = figure('Position',[624 110 900 868]);
% S = set_subplots(4,1,0.065,0.065);
% axes(S(1)) %#ok<*MAXES>
% bar(S(1),E_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
% set(S(1),'XTickLabel',num2cell(burst_types/2))
% L1 = line(xlim,[E_pulseon E_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
% L2 = line(xlim,[E_burston E_burston],'Color',[255 204 0]/255,'LineWidth',2);
% legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
% title('Efficiency')
% 
% % Plot latency and jitter
% axes(S(2))
% bar(S(2),L_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
% hold on
% errorbar(L_frequency,J_frequency,'k+')
% set(S(2),'XTickLabel',num2cell(burst_types/2))
% L1 = line(xlim,[L_pulseon L_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
% line(xlim,[L_pulseon+J_pulseon L_pulseon+J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
% line(xlim,[L_pulseon-J_pulseon L_pulseon-J_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
% L2 = line(xlim,[L_burston L_burston],'Color',[255 204 0]/255,'LineWidth',2);
% line(xlim,[L_burston+J_burston L_burston+J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
% line(xlim,[L_burston-J_burston L_burston-J_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
% legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
% title('Latency')
% 
% % Plot evoked firing rate
% axes(S(3))
% bar(S(3),M_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
% set(S(3),'XTickLabel',num2cell(burst_types/2))
% L1 = line(xlim,[M_pulseon M_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
% L2 = line(xlim,[M_burston M_burston],'Color',[255 204 0]/255,'LineWidth',2);
% L3 = line(xlim,[B_pulseon B_pulseon],'Color',[0 153 255]/255,'LineWidth',2,'LineStyle','--');
% L4 = line(xlim,[B_burston B_burston],'Color',[255 204 0]/255,'LineWidth',2,'LineStyle','--');
% legend([L1 L2 L3 L4],{'PulseOn' 'BurstOn' 'Baseline' 'Baseline'},...
%     'Location','EastOutside')
% title('Evoked firing rate')
% 
% % Plot H-index
% axes(S(4))
% bar(S(4),Hindex_frequency,'BarWidth',0.5,'EdgeColor','k','FaceColor','w')
% set(S(4),'XTickLabel',num2cell(burst_types/2))
% L1 = line(xlim,[Hindex_pulseon Hindex_pulseon],'Color',[0 153 255]/255,'LineWidth',2);
% L2 = line(xlim,[Hindex_burston Hindex_burston],'Color',[255 204 0]/255,'LineWidth',2);
% line(xlim,[0.01 0.01],'Color','r','LineWidth',2)
% legend([L1 L2],{'PulseOn' 'BurstOn'},'Location','EastOutside')
% title('SALT_Pvalue');
% filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_Efficiency_Latency_Hidx.fig'];
% saveas(gcf, filename);
% close all;

% Plot light-evoked and spont. spikes
HW = plotwaveforms(cellid,'spont',true,'evoked',true,'correlation',true,'maxnum',5000); % 'maxnum',  2000 - maximum number of spikes to plot
filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_evoked_waveform.fig'];
savefig(HW.H_evoked,filename);
filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_spont_waveform.fig'];
savefig(HW.H_spont,filename);
filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_comp_waveform.fig'];
savefig(HW.H_compare,filename);
close all;

% Plot MClust projections
HM = plot_mclust_projections2(cellid,'feature_names',{'feature_Peak','feature_Energy'},...
    'stim_period',[A1 A2]);
for i=1:length(fieldnames(HM))
filename = [resdir '\' char(regexprep(cellid,'\.','_'))  '_ClustSpace' '_' num2str(i) '.fig'];
savefig(i,filename);
end


% SDF
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, averPSTH H] = findStimPeriod(cellid);
filename = [resdir '\' char(regexprep(cellid,'\.','_'))  'SDF.fig'];
savefig(gcf,filename);
save([resdir '\TAGGEDPROP_' regexprep(cellid,'\.','_') 'AverPSTH.mat'],...
        'averPSTH');
    
close all;
% Distance from light-evoked noise
[ID_amp Lr_amp] = lightcluster(cellid,'feature_names',{'feature_Peak','feature_Energy'},...
    'stim_period',[A1 A2]);
[ID_PC Lr_PC] =lightcluster(cellid,'feature_names',{'feature_Peak','feature_WavePC1'},...
    'stim_period',[A1 A2]);
% % HD = figure;
% % axes;
% % str = {'Cluster quality measures for light-evoked spikes:';...
% %     ['ID (feature_Peak, feature_Energy): ' num2str(lID_amp)];...
% %     ['L-ratio (feature_Peak,feature_Energy): ' num2str(lLr_amp)];...
% %     ['ID (feature_WavePC1, feature_Energy): ' num2str(lID_PC)];...
% %     ['L-ratio (feature_WavePC1, feature_Energy): ' num2str(lLr_PC)];...
% %     ['preference: ' num2str(Pref)];...
% %     ['preference: ' num2str(Pref2)]};
% % uicontrol('Style','text','Unit','normalized','Position',...
% %     [0.18 0.3 0.65 0.5],'FontSize',12,'HorizontalAlignment',...
% %     'left','String',str,'BackgroundColor',[1 1 1]);
% % filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_Summary.fig'];
% % saveas(gcf, filename);
% % axis off

% BurstOn and PulseOn PSTH
HR1 = plot_raster_psth(cellid,'BurstOn',true);
filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_BurstAlign.fig'];
saveas(gcf, filename);

HR2 = plot_raster_psth(cellid,'PulseOn',true);
filename = [resdir '\' char(regexprep(cellid,'\.','_')) '_PulseAlign.fig'];
saveas(gcf, filename);
    

% Save
if issave
    save([resdir '\TAGGEDPROP_' regexprep(cellid,'\.','_') '.mat'],...
        'R',...
        'ID_amp', 'Lr_amp','ID_PC', 'Lr_PC',...
        'Hindex_act', 'D_KL_act','Hindex_inh', 'D_KL_inh',...
        'Hindex_act_burston', 'D_KL_act_burston', 'Hindex_inh_burston', 'D_KL_inh_burston',...
        'E_pulseon','L_pulseon','J_pulseon','B_pulseon','M_pulseon','lim_act1','lim_act2','Wp_act', 'MinRate', 'lim_inh1', 'lim_inh2', 'Wp_inh',... 
        'E_burston','L_burston','J_burston','B_burston','M_burston');
end
close all
