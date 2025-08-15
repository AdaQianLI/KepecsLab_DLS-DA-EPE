function figData=Online_LEDTuningPlot(action,TrialSequence,figData,thisTrial,thisNidaq)

global BpodSystem S

NbOfTrialTypes=S.NumTrialTypes;
MaxTrials=S.MaxTrials;
minxPhoto=S.GUI.TimeMin; maxxPhoto=S.GUI.TimeMax;
minyPhoto=S.GUI.NidaqMin; maxyPhoto=S.GUI.NidaqMax;
baseline=[1 20]; %Data points

switch action
    case 'ini'        
try
    close 'LED Tuning Curve';
end
%% Data initialization
for i=1:NbOfTrialTypes
	thisPower=sprintf('Power_%.0f',i);
    figData.(thisPower).X=[minxPhoto maxxPhoto];
    figData.(thisPower).Y=[0 0];
    figData.(thisPower).Data=[];
end

%% Figure initialization        
figData.figPlot=figure('Name','LED Tuning Curve','Position', [800 400 600 700], 'numbertitle','off');
ProtoSummary=sprintf('%s : %s -- %s - %s',...
    date, BpodSystem.GUIData.SubjectName, ...
    BpodSystem.GUIData.ProtocolName);
MyBox = uicontrol('style','text');
set(MyBox,'String',ProtoSummary, 'Position',[10,1,400,20]);

%% Sequence subplot
figData.Sequence.X=1:1:MaxTrials;
figData.Sequence.Y=TrialSequence;
figData.Sequence.YLabel=S.TrialsNames;
figData.Sequence.YTick=1:1:NbOfTrialTypes;

figData.Suplot(1)=subplot(3,3,[1 2],'XLim',[0 MaxTrials+1],'YLim',[0 NbOfTrialTypes+1],'YDir','reverse',...
                                    'YTick',figData.Sequence.YTick,'YTickLabel',figData.Sequence.YLabel,'YTickLabelRotation',45);
hold on;
title('Trials sequence'); xlabel('Trials'); ylabel('Power');

figData.Sequence.Plot=plot(TrialSequence,'ko');
figData.Sequence.CurTrialPlot=plot([1 1],[0 NbOfTrialTypes+1],'-r');

%% Previous Nidaq subplot
figData.Subplot(2)=subplot(3,3,3,'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]); hold on;
title('Previous recordings'); xlabel('Time(sec)'); ylabel('DF/F (%)');
figData.PreviousNidaq.Plot=plot([minxPhoto maxxPhoto],[0 0],'-k');

%% Average Subplot
figData.Subplot(3)=subplot(3,3,[4 5],'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]); hold on;
title('Average'); xlabel('Time(sec)'); ylabel('DF/F (%)');
for i=1:NbOfTrialTypes
	thisPower=sprintf('Power_%.0f',i);
	figData.(thisPower).plot=plot(figData.(thisPower).X,figData.(thisPower).Y,'-');
end
    legend(S.TrialsNames,'Location','Northeast');

%% Bleach subplot
figData.Subplot(4)=subplot(3,3,6); hold on;
title('Bleaching'); xlabel('Trial Number'); ylabel('Normalized DF/F');
figData.Bleaching.X=1:1:MaxTrials;
figData.Bleaching.Y=ones(1,MaxTrials);
figData.Bleaching.Plot=plot(figData.Bleaching.X,figData.Bleaching.Y,'og');

    case 'update'
thisTrialType=TrialSequence(thisTrial);
Time=thisNidaq(:,1);
NidaqRaw=thisNidaq(:,2);
NidaqDFF=thisNidaq(:,3);  

%% Update axes according to GUI        
set(figData.Subplot(2),'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]);
set(figData.Subplot(3),'XLim',[minxPhoto maxxPhoto],'YLim',[minyPhoto maxyPhoto]);
     
%% Update Sequence plot
set(figData.Sequence.CurTrialPlot,'XData',[thisTrial+1 thisTrial+1]);

%% Update Previous Nidaq plot
subplot(3,3,3); hold on;
title(S.TrialsNames{thisTrialType});
set(figData.PreviousNidaq.Plot,'XData',Time,'YData',NidaqDFF);
hold off

%% Updata Sounds Plot
thisPower=sprintf('Power_%.0f',thisTrialType);
if isnan(figData.(thisPower).Data)
    figData.(thisPower).Data=NidaqDFF;
else
    figData.(thisPower).Data=[figData.(thisPower).Data NidaqDFF];
end
figData.(thisPower).X=Time;
figData.(thisPower).Y=mean(figData.(thisPower).Data,2);
                   
set(figData.(thisPower).plot,'XData',figData.(thisPower).X,'YData',figData.(thisPower).Y);

%% Update Bleaching Plot
if thisTrial==1
	figData.Bleaching.Baseline=mean(NidaqRaw(baseline(1):baseline(2)));
end
figData.Bleaching.Y(thisTrial)=mean(NidaqRaw(baseline(1):baseline(2)))/figData.Bleaching.Baseline;
set(figData.Bleaching.Plot,'YData',figData.Bleaching.Y);

end
end