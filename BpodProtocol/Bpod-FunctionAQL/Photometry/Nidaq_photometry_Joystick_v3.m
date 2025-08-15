function [photometry470Data,photometry565Data,Trajectory_X,Trajectory_Y]=Nidaq_photometry_Joystick_v3(action,Param)
global nidaq S

switch action
    case 'ini'
%% NIDAQ Initialization
% Define parameters for analog inputs.
nidaq.device                   = Param.nidaqDev;
nidaq.duration                 = S.GUI.NidaqDuration;
nidaq.sample_rate              = S.GUI.NidaqSamplingRate;
nidaq.ai_channels              = {'ai0','ai1','ai2','ai3'};       % 2 channels (X-ai0 and Y-ai2 L position-ai1 R position ai3)
nidaq.ai_data                  =[];

% Define parameters for analog outputs (important for LED modulation
nidaq.ao_channels                 = {'ao0','ao1'};           % LED1 AND LED2
nidaq.ao_data                     = [];

daq.reset
daq.HardwareInfo.getInstance('DisableReferenceClockSynchronization',true); % Necessary for this Nidaq

% Set up session and channels
nidaq.session = daq.createSession('ni');


% For the photometry - Photoreceiver + TRAJECTORY
for ch = nidaq.ai_channels
    nch=addAnalogInputChannel(nidaq.session,nidaq.device,ch,'Voltage');
    nch.TerminalConfig='SingleEnded';
end

%FOR LEDs
for ch = nidaq.ao_channels
    nch=addAnalogOutputChannel(nidaq.session,nidaq.device,ch,'Voltage');
    nch.TerminalConfig='SingleEnded';
end

% Sampling rate and continuous updating (important for queue-ing ao data)
nidaq.session.Rate = nidaq.sample_rate;
nidaq.session.IsContinuous = false;
lh{1} = nidaq.session.addlistener('DataAvailable',@Nidaq_callback);
% lh{2} = nidaq.session.addlistener('DataRequired', @(src,event) src.queueOutputData(nidaq.ao_data));

    case 'WaitToStart'
%% GET NIDAQ READY TO RECORD
    nidaq.ai_data = [];
    
    nidaq.LED470              = Nidaq_modulation_joystick(S.GUI.LED470_1Amp,S.GUI.LED470_1Freq,S.GUI.LED470_1Mod); %1 channel
    if S.GUI.LED470_2Mod
    nidaq.LED565              = Nidaq_modulation_joystick(S.GUI.LED470_2Amp,S.GUI.LED470_2Freq,S.GUI.LED470_2Mod); % 2channels / 2photodetec
    elseif S.GUI.LED565_Mod
    nidaq.LED565              = Nidaq_modulation_joystick(S.GUI.LED565_Amp,S.GUI.LED565_Freq,S.GUI.LED565_Mod);
    else
    nidaq.LED565              = Nidaq_modulation_joystick(0,S.GUI.LED565_Freq,S.GUI.LED565_Mod);
    end

nidaq.ao_data           = [nidaq.LED470 nidaq.LED565];       
nidaq.session.queueOutputData(nidaq.ao_data);  
    
nidaq.session.NotifyWhenDataAvailableExceeds = nidaq.sample_rate/5;
nidaq.session.prepare();            %Saves 50ms on startup time, perhaps more for repeats.
nidaq.session.startBackground();    % takes ~0.1 second to start and release control. 

    case 'Stop'
%% STOP NIDAQ
    nidaq.session.stop()
    wait(nidaq.session) % Wait until nidaq session stop
    nidaq.session.outputSingleScan(zeros(1,length(nidaq.ao_channels))); % drop output back to 0
     case 'Save'
%% Save Data
photometry470Data=[];
photometry565Data=[];
 
% reallocates raw data

    photometry470Data = nidaq.ai_data(:,1);
    if S.GUI.DbleChanels == 1
        photometry565Data = nidaq.ai_data(:,2);
    end


    Trajectory_X      = nidaq.ai_data(:,3);
    Trajectory_Y      = nidaq.ai_data(:,4);

% saves output channels for photometry

photometry470Data = [photometry470Data nidaq.ao_data(1:size(photometry470Data,1),1)];  %col1: fluroscence col2: excitation light
if S.GUI.DbleChanels
    photometry565Data = [photometry565Data nidaq.ao_data(1:size(photometry565Data,1),2)];
else
    photometry565Data = [];
end

end     
end