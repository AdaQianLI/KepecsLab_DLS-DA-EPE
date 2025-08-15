function [photometryData,wheelData,photometry2Data]=Nidaq_photometry(action,Param)
global nidaq S

switch action
    case 'ini'
%% NIDAQ Initialization
% Define parameters for analog inputs and outputs.
nidaq.device            = Param.nidaqDev;
nidaq.duration      	= S.GUI.NidaqDuration;
nidaq.sample_rate     	= S.GUI.NidaqSamplingRate;
nidaq.ai_channels       = {'ai0','ai1'};  
nidaq.ai_data           = [];
nidaq.ao_channels       = {'ao0','ao1'};           % LED1 and LED2
nidaq.ao_data           = [];

daq.reset
daq.HardwareInfo.getInstance('DisableReferenceClockSynchronization',true); % Necessary for this Nidaq

%create nidaq session
nidaq.session = daq.createSession('ni');

% For the photometry - Photoreceiver + LEDs
for ch = nidaq.ai_channels
    nch=addAnalogInputChannel(nidaq.session,nidaq.device,ch,'Voltage');
    nch.TerminalConfig='SingleEnded';
end
for ch = nidaq.ao_channels
    nch=addAnalogOutputChannel(nidaq.session,nidaq.device,ch,'Voltage');
    nch.TerminalConfig='SingleEnded';
end

% For the running wheel - quadrature rotary encoder
% nidaq.counter = addCounterInputChannel(nidaq.session,nidaq.device,0,'Position');
% nidaq.counter.EncoderType = 'X1';
% nidaq.counter.ZResetEnable = true;
% nidaq.counter.ZResetCondition = 'BothLow';
% nidaq.counter.ZResetValue = 0;

% Sampling rate
nidaq.session.Rate = nidaq.sample_rate;
nidaq.session.IsContinuous = false;
lh{1} = nidaq.session.addlistener('DataAvailable',@Nidaq_callback);

    case 'WaitToStart'
%% GET NIDAQ READY TO RECORD
nidaq.ai_data            = [];
if S.GUI.Photometry
    nidaq.LED1              = Nidaq_modulation(S.GUI.LED1_Amp,S.GUI.LED1_Freq); %1 channel
    if S.GUI.DbleFibers == 0
    nidaq.LED2              = Nidaq_modulation(S.GUI.LED2_Amp,S.GUI.LED2_Freq); %2 channel / 1photodetec/565
        else
    nidaq.LED2              = Nidaq_modulation(S.GUI.LED1b_Amp,S.GUI.LED1b_Freq); % 2channels / 2photodetec
    end
else
    nidaq.LED1              = Nidaq_modulation(0,S.GUI.LED1_Freq);
    nidaq.LED2              = Nidaq_modulation(0,S.GUI.LED2_Freq);
end
nidaq.ao_data           = [nidaq.LED1 nidaq.LED2];
nidaq.session.queueOutputData(nidaq.ao_data);

nidaq.session.NotifyWhenDataAvailableExceeds = nidaq.sample_rate/5;
nidaq.session.prepare();
nidaq.session.startBackground();

    case 'Stop'
%% STOP NIDAQ
    nidaq.session.stop()
    wait(nidaq.session) % Wait until nidaq session stop
    nidaq.session.outputSingleScan(zeros(1,length(nidaq.ao_channels))); % drop output back to 0 
    case 'Save'
%% Save Data
photometryData=[];
photometry2Data=[];
wheelData=[];

% reallocates raw data
if S.GUI.Photometry
    photometryData = nidaq.ai_data(:,1);
    if  S.GUI.DbleFibers || S.GUI.RedChannel 
        photometry2Data = nidaq.ai_data(:,2);
    end
end
if S.GUI.Wheel
    wheelData      = nidaq.ai_data(:,3);
end

% saves output channels for photometry
% if S.GUI.Photometry
%     if S.GUI.Modulation
%         if S.GUI.LED2_Amp==0 % LED2_Amp channel / 1photodetec/565
%             photometryData = [photometryData nidaq.ao_data(1:size(photometryData,1),1)];
%             if S.GUI.DbleFibers == 1
%                 photometry2Data = [photometry2Data nidaq.ao_data(1:size(photometryData,1),2)];
%             else
%                 photometry2Data = [];
%             end
%         else
%             photometryData = [photometryData nidaq.ao_data(1:size(photometryData,1),:)];
%         end
%     end
% end  
end     
end