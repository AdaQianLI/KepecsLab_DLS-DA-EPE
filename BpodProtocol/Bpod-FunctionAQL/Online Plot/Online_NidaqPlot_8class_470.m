function figData=Online_NidaqPlot_8class_470(action,Phase,trialsResNames,figData,newData470,nidaqRaw,curTrialType)
global BpodSystem S
%% general ploting parameters

labelx='Time (sec)'; 
labely='DF/F'; 
minx=S.GUI.TimeMin_photometry; 
maxx=S.GUI.TimeMax_photometry;
xstep=1;   
xtickvalues=minx:xstep:maxx;
miny=S.GUI.NidaqMin; 
maxy=S.GUI.NidaqMax;
MeanThickness=2;

switch action
    case 'ini'
%% Close pre-existing plot and test parameters
try
    close 'Online Nidaq Plot_470nm';
end
%% Create Figure
figPlot=figure('Name','Online Nidaq Plot_470nm','Position', [550 50 600 700], 'numbertitle','off'); %[left bottom width height]
hold on
ProtoSummary=sprintf('%s : %s -- %s - %s',...
    date, BpodSystem.GUIData.SubjectName, ...
    BpodSystem.GUIData.ProtocolName, Phase);
ProtoLegend=uicontrol('style','text');
set(ProtoLegend,'String',ProtoSummary); 
set(ProtoLegend,'Position',[10,1,400,20]);

trialsNames_8class=trialsResNames;

%% Current Nidaq plot
% lastsubplot=subplot(3,2,[1 2]);
lastsubplot=subplot(5,2,[1 2]);
hold on
title('Nidaq 470 recording');
xlabel(labelx); 
ylabel('Voltage');
ylim auto;
set(lastsubplot,'XLim',[minx maxx],'XTick',xtickvalues);%'YLim',[miny maxy]
lastplotRaw=plot([-5 5],[0 0],'-k');
lastplot470=plot([-5 5],[0 0],'-g','LineWidth',MeanThickness);
%lastplot565=plot([-5 5],[0 0],'-r','LineWidth',MeanThickness);
hold off

%% Plot previous recordings
for i=1:size(trialsNames_8class,2)
%     photosubplot(i)=subplot(3,2,i+2);
    photosubplot(i)=subplot(5,2,i+2);
    hold on
    xlabel(labelx); 
    ylabel(labely);
    ylim auto;
    set(photosubplot,'XLim',[minx maxx],'XTick',xtickvalues);
    rewplot(i)=plot([0 0],[-1,1],'-b');
    meanplot(i)=plot([-5 5],[0 0],'-r');
    title(trialsNames_8class{1,i});
    hold off

set(photosubplot(i),'XLabel',[],'YLabel',[]);
end

%Save the figure properties
figData.fig=figPlot;
figData.lastsubplot=lastsubplot;
figData.lastplotRaw=lastplotRaw;
figData.lastplot470=lastplot470;
%figData.lastplot565=lastplot565;
figData.photosubplot=photosubplot;
figData.meanplot=meanplot;

    case 'update'
        
set(figData.lastplotRaw, 'Xdata',nidaqRaw(:,1),'YData',nidaqRaw(:,2));
set(figData.lastplot470, 'Xdata',newData470(:,1),'YData',newData470(:,2));
%set(figData.lastplot565, 'Xdata',newData565(:,1),'YData',newData565(:,2));        
        
         if curTrialType<=8
%% Compute new average trace
allData=get(figData.photosubplot(curTrialType), 'UserData');
dataSize=size(allData,2);
allData(:,dataSize+1)=newData470(:,3);
set(figData.photosubplot(curTrialType), 'UserData', allData);
meanData=mean(allData,2);
factor=50;
curSubplot=figData.photosubplot(curTrialType);
set(figData.meanplot(curTrialType), 'XData',newData470(:,1),'YData',meanData,'LineWidth',MeanThickness);
set(curSubplot,'NextPlot','add');
plot(newData470(1:factor:end,1),newData470(1:factor:end,3),'-k','parent',curSubplot);
uistack(figData.meanplot(curTrialType), 'top');
hold off
         end
%% Update GUI plot parameters
 for i=1:size(trialsResNames,2)
     set(figData.photosubplot(i),'XLim',[minx maxx],'XTick',xtickvalues)
end
end
end