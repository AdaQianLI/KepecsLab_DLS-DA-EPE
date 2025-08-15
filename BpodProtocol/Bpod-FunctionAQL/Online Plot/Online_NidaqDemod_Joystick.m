function [NidaqDemod, NidaqRaw]=Online_NidaqDemod_Joystick(rawData,refData,modFreq,modAmp,StateToZero,thisTrial,mod)
global BpodSystem S

decimateFactor=S.GUI.DecimateFactor_photon;
duration=S.GUI.NidaqDuration;
sampleRate=S.GUI.NidaqSamplingRate;
baseline_begin=S.GUI.BaselineBegin;
baseline_end=S.GUI.BaselineEnd;
lowCutoff=20;
pad=1;



%% Filter
    lowCutoff = lowCutoff/(sampleRate/2); % normalized CutOff by half SampRate (see doc)
    [b, a,k] = butter(10, lowCutoff, 'low'); 
    [sos, g] = zp2sos(b,a,k);
 if mod  
     %% Prepare reference data and generates 90deg shifted ref data
refData             = refData(1:length(rawData),1);   % adjust length of refData to rawData
refData             = refData-mean(refData);          % suppress DC offset
samplesPerPeriod    = (1/modFreq)/(1/sampleRate);
quarterPeriod       = round(samplesPerPeriod/4);
refData90           = circshift(refData,[1 quarterPeriod]);

%% Quadrature decoding and filtering
processedData_0     = rawData .* refData;
processedData_90    = rawData .* refData90;
    if pad == 1 % pad the data to suppress windows effect upon filtering
        paddedData_0        = processedData_0(1:sampleRate, 1);
        paddedData_90       = processedData_0(1:sampleRate, 1);
        demodDataFilt_0     = filtfilt(sos,g,[paddedData_0; processedData_0]);
        demodDataFilt_90    = filtfilt(sos,g,[paddedData_90; processedData_90]);        
        processedData_0     = demodDataFilt_0(sampleRate + 1: end, 1);
        processedData_90    = demodDataFilt_90(sampleRate + 1: end, 1);
    else
        processedData_0     = filtfilt(b,a,processedData_0);
        processedData_90    = filtfilt(b,a,processedData_90); 
    end
    
demodData = (processedData_0 .^2 + processedData_90 .^2) .^(1/2);

%% Correct for amplitude of reference
demodData=demodData*2/modAmp;
else
%     demodData=rawData;
            paddedData = rawData(1:min(sampleRate, size(rawData, 1)), 1);
            temp = filtfilt(sos, g, [paddedData; rawData]);
            demodData = temp(length(paddedData) + 1:end, 1);
            demodData=demodData*2;
end
%% Expeced Data set
SampRate=sampleRate/decimateFactor;
ExpectedSize=duration*SampRate;
Data=NaN(ExpectedSize,1);
TempData=decimate(demodData,decimateFactor);
Data(1:length(TempData))=TempData;

%% DF/F calculation
Fbaseline=mean(Data(baseline_begin*SampRate:baseline_end*SampRate));
DFF=100*(Data-Fbaseline)/Fbaseline;

%% Time
Time=linspace(0,duration,ExpectedSize);
TimeToZero=BpodSystem.Data.RawEvents.Trial{1,thisTrial}.States.(StateToZero)(1,1);
Time=Time'-TimeToZero;

%% Raw Data
ExpectedSizeRaw=duration*sampleRate;
DataRaw=NaN(ExpectedSizeRaw,1);
DataRaw(1:length(rawData))=rawData;

TimeRaw=linspace(0,duration,ExpectedSizeRaw);
TimeRaw=TimeRaw'-TimeToZero;
%% NewDataSet
NidaqDemod(:,1)=Time;
NidaqDemod(:,2)=Data;
NidaqDemod(:,3)=DFF;

NidaqRaw(:,1)=TimeRaw;
NidaqRaw(:,2)=DataRaw;
end