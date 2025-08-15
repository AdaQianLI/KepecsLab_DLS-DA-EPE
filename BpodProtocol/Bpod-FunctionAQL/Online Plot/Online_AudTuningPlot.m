function figData=Online_AudTuningPlot(action,TrialSequence,figData,thisTrial,thisNidaq)

global BpodSystem S

NbOfTrialTypes=S.NumTrialTypes;
MaxTrials=S.MaxTrials;
minxPhoto=S.GUI.TimeMin; maxxPhoto=S.GUI.TimeMax;
minyPhoto=S.GUI.NidaqMin; maxyPhoto=S.GUI.NidaqMax;
baseline=[1 20]; %Data points
soundResponse=[S.GUI.CueBegin S.GUI.CueEnd];

switch action
    case 'ini'        
try
    close 'Auditory Tuning Curve';
end
%% Data initialization
for i=1:NbOfTrialTypes
	thisSound=sprintf('Sound_%.0d',i);
    figData.(thisSound).X=[minxPhoto maxxPhoto];
    figData.(thisSound).Y=[0 0];
    figData.(thisSound).Data=[];
    figData.Tuning.X(i)=i;
    figData.Tuning.Y(i)=NaN;
end

%% Figure initialization        
figData.figPlot=figure('Name','Auditory Tuning Curve','Position', [800 400 600 700], 'numbertitle','off');
ProtoSummary=sprintf('%s : %s -- %s - %s',...
    date, BpodSystem.GUIData.SubjectName, ...
    BpodSystem.GUIData.ProtocolName);
MyBox = uicontrol('style','text');
set(MyBox,'String',ProtoSummary, 'Position',[10,1,400,20]);

%% Sequence subplot
figData.Sequence.X=1:1:MaxTrials;
figData.Sequence.Y=TrialSequence;

figData.Sequence.Ylabel={'White Noise','Sweep','Pure Tones'};
figData.Sequence.YTick=[1 2.5 4.5];

figData.Suplot(1)=subplot(3,3,[1 2],'XLim',[0 MaxTrials+1],'YLim',[0 NbOfTrialTypes+1],'YDir','reverse',...
                                    'YTick',figData.Sequence.YTick,'YTickLabel',figData.Sequence.Ylabel,'YTickLabelRotation', 45 );
hold on;
title('Trials sequence'); xlabel('Trials'); ylabel('Sounds');

figData.Sequence.Plot=plot(TrialSequence,'ko');
figData.Sequence.PlotLine1=plot([1 MaxTrials],[1.5 1.5],'-g');
figData.Sequence.PlotLine2=plot([1 MaxTrials],[3.5 3.5],'-g');
figData.Sequence.CurTrialPlot=plot([1 1],[0 NbOfTrialTypes+1],'-r');

%% Previous Nidaq subplot
figData.Subplot(2)=subplot(3,3,3,'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]); hold on;
title('Previous recordings'); xlabel('Time(sec)'); ylabel('DF/F (%)');
figData.PreviousNidaq.Plot=plot([minxPhoto maxxPhoto],[0 0],'-k');

%% WhiteNoise subplot
figData.Subplot(3)=subplot(3,3,4,'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]); hold on;
title('WhiteNoise'); xlabel('Time(sec)'); ylabel('DF/F (%)');
if S.GUI.WhiteNoise
    figData.Sound_1.plot=plot(figData.Sound_1.X,figData.Sound_1.Y,'-k');
    legend(S.TrialsNames{1},'Location','Northeast');
end

%% Sweeps subplot
figData.Subplot(4)=subplot(3,3,5,'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]); hold on;
title('Sweeps'); xlabel('Time(sec)'); ylabel('DF/F (%)');
if S.GUI.Sweeps
    figData.Sound_2.plot=plot(figData.Sound_2.X,figData.Sound_2.Y,'-r');
	figData.Sound_3.plot=plot(figData.Sound_3.X,figData.Sound_3.Y,'-b');
    legend(S.TrialsNames{2},S.TrialsNames{3},'Location','Northeast'); 
end

%% Pure Tones subplot
figData.Subplot(5)=subplot(3,3,6,'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]); hold on;
title('Pure Tones'); xlabel('Time(sec)'); ylabel('DF/F (%)');
if S.GUI.PureTones
    for i=4:NbOfTrialTypes
        thisSound=sprintf('Sound_%.0d',i);
        figData.(thisSound).plot=plot(figData.(thisSound).X,figData.(thisSound).Y,'-');
    end
    legend(S.TrialsNames{4:NbOfTrialTypes},'Location','Northeast');
end

%% Bleach subplot
figData.Subplot(6)=subplot(3,3,7); hold on;
title('Bleaching'); xlabel('Trial Number'); ylabel('Normalized DF/F');
figData.Bleaching.X=1:1:MaxTrials;
figData.Bleaching.Y=ones(1,MaxTrials);
figData.Bleaching.Plot=plot(figData.Bleaching.X,figData.Bleaching.Y,'og');

%% Auditory Tuning subplot
figData.Subplot(7)=subplot(3,3,[8 9],'XLim',[0 NbOfTrialTypes+1],'YLim',[minyPhoto maxyPhoto],...
                                        'XTick',figData.Tuning.X,'XTickLabel',S.TrialsNames,'XTickLabelRotation', 45); hold on;
title('Auditory Tuning'); ylabel('DF/F (%)');
figData.Tuning.Plot=plot(figData.Tuning.X,figData.Tuning.Y,'sb');
figData.Tuning.PlotLine1=plot([1.5 1.5],[minyPhoto maxyPhoto],'-b');
figData.Tuning.PlotLine2=plot([3.5 3.5],[minyPhoto maxyPhoto],'-b');

    case 'update'
thisTrialType=TrialSequence(thisTrial);
Time=thisNidaq(:,1);
NidaqRaw=thisNidaq(:,2);
NidaqDFF=thisNidaq(:,3);  

%% Update axes according to GUI        
for i=2:5
    set(figData.Subplot(i),'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]);
end
set(figData.Subplot(7),'YLim',[minyPhoto maxyPhoto]);
        
%% Update Sequence plot
set(figData.Sequence.CurTrialPlot,'XData',[thisTrial+1 thisTrial+1]);

%% Update Previous Nidaq plot
subplot(3,3,3); hold on;
title(S.TrialsNames{thisTrialType});
set(figData.PreviousNidaq.Plot,'XData',Time,'YData',NidaqDFF);
hold off

%% Updata Sounds Plot
thisSound=sprintf('Sound_%.0d',thisTrialType);
if isnan(figData.(thisSound).Data)
    figData.(thisSound).Data=NidaqDFF;
else
    figData.(thisSound).Data=[figData.(thisSound).Data NidaqDFF];
end
figData.(thisSound).X=Time;
figData.(thisSound).Y=mean(figData.(thisSound).Data,2);
                   
set(figData.(thisSound).plot,'XData',figData.(thisSound).X,'YData',figData.(thisSound).Y);

%% Update Bleaching Plot
if thisTrial==1
	figData.Bleaching.Baseline=mean(NidaqRaw(baseline(1):baseline(2)));
end
figData.Bleaching.Y(thisTrial)=mean(NidaqRaw(baseline(1):baseline(2)))/figData.Bleaching.Baseline;
set(figData.Bleaching.Plot,'YData',figData.Bleaching.Y);

%% Update Tuning curve Plot
for i=1:NbOfTrialTypes
    thisSound=sprintf('Sound_%.0d',i);
    if ~isempty(figData.(thisSound).Data)
        figData.Tuning.Y(i)=mean(figData.(thisSound).Y(Time>soundResponse(1) & Time<soundResponse(2)));
    end
end
set(figData.Tuning.Plot,'YData',figData.Tuning.Y);
end

end