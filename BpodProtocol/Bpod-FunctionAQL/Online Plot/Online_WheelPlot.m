function figData=Online_WheelPlot(action,figData,DataCount,StateToZero,thisTrial,thisLicks)
global BpodSystem S
%% Plot Parameters
labelx='Time(sec)'; labely='AngularSpeed deg/sec'; labelyy='Angle deg';
minx=S.GUI.TimeMin; maxx=S.GUI.TimeMax; xstep=1; xtickvalues=minx:xstep:maxx;

switch action
    case 'ini'
try
    close 'Online Wheel';
end
% Figure
ScrSze=get(0,'ScreenSize');
FigSze=[ScrSze(3)*2/3 ScrSze(2)+40 ScrSze(3)*1/3 300];
figPlot=figure('Name','Online Wheel','Position',FigSze,'numbertitle','off');
hold on
ProtoSummary=sprintf('%s : %s -- %s - %s',...
    date, BpodSystem.GUIData.SubjectName, ...
    BpodSystem.GUIData.ProtocolName);
ProtoLegend=uicontrol('style','text');
set(ProtoLegend,'String',ProtoSummary); 
set(ProtoLegend,'Position',[10,1,400,20]);

% Plot
wheelSubplot=subplot(1,1,1);
hold on
%yyaxis left
SpeedPlot=plot(0,0,'-k');
xlabel(labelx);ylabel(labely);

%yyaxis right
AnglePlot=plot(0,0,'-r');
%ylabel(labelyy);
LickPlot=scatter(0,0,10,'v','filled');
hold off

% Save the figure handle
figData.fig=figPlot;
figData.subplot=wheelSubplot;
figData.SpeedPlot=SpeedPlot;
figData.AnglePlot=AnglePlot;
figData.LickPlot=LickPlot;

    case 'update'
%% Parameters
SamplingRate    = S.GUI.NidaqSamplingRate;
DecimateFactor  = 61;
CounterNBits    = 32;
EncoderCPR      = 1024;
newSR           = SamplingRate/DecimateFactor; % 10 msec for decimate=61
deltaT          = 1/newSR;

%% Decoding
signedThreshold = 2^(CounterNBits-1);
signedData = DataCount;
signedData(signedData > signedThreshold) = signedData(signedData > signedThreshold) - 2^CounterNBits;
DataDeg  = signedData * 360/EncoderCPR;
DataDeg  = decimate((DataDeg),DecimateFactor);
ActualDuration=length(DataDeg)*deltaT-deltaT;

%% Speed
%time array
Time=deltaT:deltaT:ActualDuration;
TimeToZero=BpodSystem.Data.RawEvents.Trial{1,thisTrial}.States.(StateToZero)(1,1);
Time=Time'-TimeToZero;

%speed calculation
% speed=diff(DataDeg)/deltaT;

%% Update Plot
%hold on
%set(figData.SpeedPlot,'Xdata',Time,'YData',speed);
set(figData.AnglePlot,'Xdata',Time,'YData',DataDeg(1:length(Time)));
ylim auto
LimY=get(figData.subplot,'YLim');
LicksY=ones(length(thisLicks),1)*LimY(2);
set(figData.LickPlot,'Xdata',thisLicks,'YData',LicksY);
set(figData.subplot,'XLim',[minx maxx],'XTick',xtickvalues);

end
end