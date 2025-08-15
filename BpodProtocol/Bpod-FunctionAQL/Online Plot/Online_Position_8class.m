function figData=Online_Position_8class(action,Phase,trialsResNames,figData,newData,curTrialType)
global BpodSystem S
%% general ploting parameters
labelx='X(v)'; 
labely='Y(v)'; 
MeanThickness=2;

switch action
    case 'ini'
%% Close pre-existing plot and test parameters
try
    close 'Online Position Plot';
end
%% Create Figure
figPlot=figure('Name','Online Position Plot','Position', [10 50 350 1000], 'numbertitle','off');  %[left bottom width height]
hold on
ProtoSummary=sprintf('%s : %s -- %s - %s',...
    date, BpodSystem.GUIData.SubjectName, ...
    BpodSystem.GUIData.ProtocolName, Phase);
ProtoLegend=uicontrol('style','text');
set(ProtoLegend,'String',ProtoSummary); 
set(ProtoLegend,'Position',[10,1,400,20]);

%% Current Nidaq plot
lastsubplot=subplot(5,2,1);
hold on
title('Position recording');
xlabel(labelx); 
ylabel(labely);
all_x=[];
all_y=[];
xlim([0 5]);
ylim([0 5]);
hold off

trialsNames_8class=trialsResNames;


%% Plot previous recordings
for i=1:size(trialsNames_8class,2)
    motorsubplot(i)=subplot(5,2,i+2);
    hold on
    xlabel(labelx); 
    ylabel(labely);
    set(motorsubplot(i),'XLim',[0 5],'YLim',[0 5]);
    meanplot(i)=plot([0 0],[0 0],'-r'); 
    motorplot(i)=plot([0 0],[-1,1],'-b');
    set(motorplot, 'XData',[],'YData',[]);
    title(trialsNames_8class{1,i});
    hold off
    set(motorsubplot(i),'XLabel',[],'YLabel',[]);
end


%Save the figure properties
figData.fig=figPlot;
figData.lastsubplot=lastsubplot;
figData.motorsubplot=motorsubplot;
figData.meanplot=meanplot;
figData.motorplot=motorplot;

    case 'update'

curSubplot0=figData.lastsubplot;
plot(newData(:,2),newData(:,3),'-g','parent',curSubplot0);
set(curSubplot0,'XLim',[0 5],'YLim',[0 5]);%'YLim',[miny maxy]

%% Compute new average trace
 all_x=get(figData.motorplot(curTrialType), 'XData')';
 all_y=get(figData.motorplot(curTrialType), 'YData')';
 all_x=[all_x newData(:,2)];
 all_y=[all_y newData(:,3)];

meanData_x=mean(all_x,2);
meanData_y=mean(all_y,2);


curSubplot=figData.motorsubplot(curTrialType);
set(curSubplot,'NextPlot','add');
plot(newData(:,2),newData(:,3),'-k','parent',curSubplot);

set(figData.motorplot(curTrialType),'XData',newData(:,2),'YData',newData(:,3));  %% handle for trace

set(figData.meanplot(curTrialType),'XData',meanData_x,'YData',meanData_y,'LineWidth',MeanThickness);

% set(figData.motorsubplot(curTrialType),'XLim',[xmin xmax],'YLim',[ymin ymax]);

uistack(figData.meanplot(curTrialType), 'top');
hold off
% %% Update GUI plot parameters
%  for i=1:8
%      
% end
end
end