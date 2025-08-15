function NidaqData=Online_NidaqPosition(Traj_X,Traj_Y,thisTrial,StateToZero)
global BpodSystem  S

%Parameters
DecimateFactor=S.GUI.DecimateFactor_traj;
Duration=S.GUI.NidaqDuration;
SampRate=S.GUI.NidaqSamplingRate;
% BaselineBegin=S.GUI.BaselineBegin;
% BaselineEnd=S.GUI.BaselineEnd;
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
% Fbaseline=mean(DataX(BaselineBegin*SampRate:BaselineEnd*SampRate));
% DFF=100*(DataX-Fbaseline)/Fbaseline;

%Time

Time=linspace(0,Duration,ExpectedSize);
TimeToZero=BpodSystem.Data.RawEvents.Trial{1,thisTrial}.States.(StateToZero)(1,1);
% TimeToEnd=BpodSystem.Data.RawEvents.Trial{1,thisTrial}.States.(StateToEnd)(1,1);
    Time=Time'-TimeToZero;
%NewDataXSet
datalenMin=round(TimeToZero*SampRate)+1;
datalenMax=round((TimeToZero+5)*SampRate);
NidaqData(:,1)=Time(datalenMin:datalenMax-1);
NidaqData(:,2)=DataX(datalenMin:datalenMax-1);
NidaqData(:,3)=DataY(datalenMin:datalenMax-1);
  
end