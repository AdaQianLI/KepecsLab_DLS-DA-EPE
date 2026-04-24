function NidaqData=Online_NidaqPosition(Traj_X,Traj_Y,thisTrial,StateToZero)
% ONLINE_NIDAQPOSITION
% Downsample and align joystick X/Y signals to the chosen task event for
% online trajectory plotting.
global BpodSystem  S

% Parameters for downsampling and alignment
DecimateFactor=S.GUI.DecimateFactor_traj;
Duration=S.GUI.NidaqDuration;
SampRate=S.GUI.NidaqSamplingRate;
SampRate=SampRate/DecimateFactor; %Hz

%Generates DataX
ExpectedSize=Duration*SampRate;
DataX=NaN(ExpectedSize,1);
DataY=NaN(ExpectedSize,1);
TempDataX=decimate(Traj_X,DecimateFactor);
TempDataY=decimate(Traj_Y,DecimateFactor);

%DataX processing
DataX(1:length(TempDataX))=TempDataX;
DataY(1:length(TempDataX))=TempDataY;

% Time axis aligned to the selected event

Time=linspace(0,Duration,ExpectedSize);
TimeToZero=BpodSystem.Data.RawEvents.Trial{1,thisTrial}.States.(StateToZero)(1,1);
    Time=Time'-TimeToZero;
%NewDataXSet
datalenMin=round(TimeToZero*SampRate)+1;
datalenMax=round((TimeToZero+5)*SampRate);
NidaqData(:,1)=Time(datalenMin:datalenMax-1);
NidaqData(:,2)=DataX(datalenMin:datalenMax-1);
NidaqData(:,3)=DataY(datalenMin:datalenMax-1);
  
end