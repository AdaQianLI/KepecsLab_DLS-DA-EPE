function figData=Online_NidaqPlot(action,subPlotTitles,Phase,figData,newData470,newData405,nidaqRaw,curTrialType)
global BpodSystem S
%% general ploting parameters
labelx='Time (sec)'; labely='DF/F'; 
minx=S.GUI.TimeMin; maxx=S.GUI.TimeMax;  xstep=1;    xtickvalues=minx:xstep:maxx;
miny=S.GUI.NidaqMin; maxy=S.GUI.NidaqMax;
MeanThickness=2;

switch action
    case 'ini'
%% Close pre-existing plot and test parameters
if S.GUI.DbleFibers == 0
    try
        close 'Online Nidaq Plot';
    end
end
%% Create Figure
ScrSze=get(0,'ScreenSize');
FigSze=[ScrSze(3)*1/3 ScrSze(2)+40 ScrSze(3)*1/3 ScrSze(4)-120];
figPlot=figure('Name','Online Nidaq Plot','Position',FigSze, 'numbertitle','off');
hold on
ProtoSummary=sprintf('%s : %s -- %s - %s',...
    date, BpodSystem.GUIData.SubjectName, ...
    BpodSystem.GUIData.ProtocolName, Phase);
ProtoLegend=uicontrol('style','text');
set(ProtoLegend,'String',ProtoSummary); 
set(ProtoLegend,'Position',[10,1,400,20]);

%% Current Nidaq plot
lastsubplot=subplot(4,2,[1 2]);
hold on
title('Nidaq recording');
xlabel(labelx); ylabel('Voltage');
ylim auto;
set(lastsubplot,'XLim',[minx maxx],'XTick',xtickvalues);%'YLim',[miny maxy]
lastplotRaw=plot([-5 5],[0 0],'-k');
lastplot470=plot([-5 5],[0 0],'-g','LineWidth',MeanThickness);
lastplot405=plot([-5 5],[0 0],'-r','LineWidth',MeanThickness);
hold off

%% Plot previous recordings
for i=1:6
    photosubplot(i)=subplot(4,2,i+2);
    hold on
    title(subPlotTitles(i));
    xlabel(labelx); ylabel(labely);
    ylim auto;
    set(photosubplot(i),'XLim',[minx maxx],'XTick',xtickvalues,'YLim',[miny maxy]);
    rewplot(i)=plot([0 0],[-1,1],'-b');
    meanplot(i)=plot([-5 5],[0 0],'-r');
    hold off
end

set(photosubplot(1),'XLabel',[]);
set(photosubplot(2),'XLabel',[],'YLabel',[]);
set(photosubplot(3),'XLabel',[]);
set(photosubplot(4),'XLabel',[],'YLabel',[])
%set(photosubplot(5),'XLabel',labelx,'YLabel',labely);
set(photosubplot(6),'YLabel',[]);

%Save the figure properties
figData.fig=figPlot;
figData.lastsubplot=lastsubplot;
figData.lastplotRaw=lastplotRaw;
figData.lastplot470=lastplot470;
figData.lastplot405=lastplot405;
figData.photosubplot=photosubplot;
figData.meanplot=meanplot;

    case 'update'
        %% Update last recording plot
set(figData.lastplotRaw, 'Xdata',nidaqRaw(:,1),'YData',nidaqRaw(:,2));
set(figData.lastplot470, 'Xdata',newData470(:,1),'YData',newData470(:,2));
set(figData.lastplot405, 'Xdata',newData405(:,1),'YData',newData405(:,2));
         if curTrialType<=6
%% Compute new average trace
allData=get(figData.photosubplot(curTrialType), 'UserData');
dataSize=size(allData,2);
allData(:,dataSize+1)=newData470(:,3);
set(figData.photosubplot(curTrialType), 'UserData', allData);
meanData=mean(allData,2);

curSubplot=figData.photosubplot(curTrialType);
set(figData.meanplot(curTrialType), 'Xdata',newData470(:,1),'YData',meanData,'LineWidth',MeanThickness);
set(curSubplot,'NextPlot','add');
plot(newData470(:,1),newData470(:,3),'-k','parent',curSubplot);
uistack(figData.meanplot(curTrialType), 'top');
hold off
         end
%% Update GUI plot parameters
 for i=1:6
     set(figData.photosubplot(i),'XLim',[minx maxx],'XTick',xtickvalues,'YLim',[miny maxy])
end
end
end