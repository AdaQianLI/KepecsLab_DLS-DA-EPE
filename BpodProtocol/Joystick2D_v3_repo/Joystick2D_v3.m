function Joystick2D_v3_noRwdcueNew %% Joystick2D_v3_repo

% JOYSTICK2D_V3
% Main Bpod protocol for joystick reaching with optional force fields,
% photometry acquisition, and online extraction of dF/F & plotting.
%
% Trial structure:
%   4 s pre-cue period -> 0.5 s cue -> 1 s delay -> 5.5 s response window
% Successful trials then enter a 1 s pre-reward delay
% Premature reaches before the go period trigger the warning
% buzzer and terminate the trial.
% -> 3-8 s ITI

global BpodSystem nidaq S

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load any previously saved GUI settings for this protocol.
ParamPC=BpodParam_PCdep_Joystick_v3();
if isempty(fieldnames(S))  % First launch: populate the GUI structure with protocol defaults.
    Bpod_Joystick2D_v3_parameter(ParamPC);
end

% Initialize the parameter GUI and pause once so the user can inspect settings before starting.
BpodParameterGUI('init', S);
BpodSystem.Pause=1;
HandlePauseCondition;
 S = BpodParameterGUI('sync', S);
 
%% Build task sounds and register them with PsychToolbox
% Sound IDs used later by SoftCodeHandler_PlaySound:

%   1 = early-reach warning buzzer (three brief beeps)
%   2 = white-noise timeout/fail sound

TooEealyCue=beepsound(S.GUI.SF,S.GUI.BeepFre,0.1);
TimeOutCue = WhiteNoise(S.GUI.SF,S.GUI.SoundDur);  %noise

% Program sound server
PsychToolboxSoundServer('init');
PsychToolboxSoundServer('Load', 1, TooEealyCue); 
PsychToolboxSoundServer('Load', 2, TimeOutCue);  %noise 
BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler_PlaySound';


 [S.TrialsMatrix, S.TrialsNames, S.TrialsResNames]=Joystick2D_v3_Phase(S.Names.Phase{S.GUI.Phase});  % Load the selected training/testing phase.
 tempSequence=WeightedRandomTrials(S.TrialsMatrix(:,2)', S.GUI.MaxTrials); % Weighted draw based on phase-specific probabilities.


 S.TrialSequence=[tempSequence]; % Planned sequence of trial identities for the session.
 
% Map each trial identity onto the force threshold value sent over Serial1.
for ii=1:length(S.TrialSequence)
    switch S.TrialSequence(ii)
        case 1
            S.force.thr(ii)=S.GUI.L0Force;
        case 2
            S.force.thr(ii)=S.GUI.L1Force;
        case 3
            S.force.thr(ii)=S.GUI.L2Force;  
        case 4
            S.force.thr(ii)=S.GUI.L3Force; 
    end
end
if S.GUI.OptoStim == 1
    PulsePal('COM5');
    load('C:\Users\Kepecs\Documents\Data\Ada\Bpod\Protocols\Joystick2D_v3_repo\LightTrain_500ms.mat');
    ProgramPulsePal(ParameterMatrix);
end
%% Sample the randomized ITI component
% The code uses a truncated exponential between 1 and 6 s. Together with the
% fixed 2 s Break state, this yields a total between-trial interval of 3-8 s.
n=length(S.TrialSequence);
lambda=1/3;
a=1;
b=6;
fa=1-exp(-lambda*a);
fb=1-exp(-lambda*b);
u=fa+(fb-fa)*rand(n,1);
iti=-log(1-u)/lambda;
iti=iti';

%% Initialize NI-DAQ acquisition and online figures

Nidaq_photometry_Joystick_v3('ini',ParamPC);
FigLick=Online_LickinJoystick('ini',S.TrialSequence,S.TrialsMatrix,S.Names.Phase{S.GUI.Phase},S.TrialsNames,S.TrialsResNames);
FigNidaq470=Online_NidaqPlot_8class_470('ini',S.Names.Phase{S.GUI.Phase},S.TrialsResNames); %'Air_Stim_Pairing'
if S.GUI.DbleChanels == 1
FigNidaq565=Online_NidaqPlot_8class_565('ini',S.Names.Phase{S.GUI.Phase},S.TrialsResNames); %'Air_Stim_Pairing'
end    
FigTraj=Online_Position_8class('ini',S.Names.Phase{S.GUI.Phase},S.TrialsResNames);


%% Main trial loop
for currentTrial = 1:length(S.TrialSequence)
S = BpodParameterGUI('sync', S); % Pull in any user edits made in the GUI before this trial starts.
      
%% Assemble the Bpod state matrix for the current trial
 	sma = NewStateMatrix();
    sma = SetGlobalTimer(sma,1,5.5); % 5.5 s response window after the go signal.

    sma = AddState(sma, 'Name','Dummy1',...
        'Timer',4,...
        'StateChangeConditions',{'Tup','CueDeliver'},...
        'OutputActions',{}); % 4 s pre-cue period described in Methods.
    sma = AddState(sma, 'Name','CueDeliver',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','Delay0','Wire3High','TooEarly'},...
        'OutputActions',{'BNCState',S.TrialsMatrix(S.TrialSequence(currentTrial),9)});   % Visual cue / external cue trigger.
    sma = AddState(sma, 'Name','Delay0',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','Delay','Wire3High','TooEarly'},... %,'Wire3High','TooEarly'
        'OutputActions',{});
    sma = AddState(sma, 'Name','Delay',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','Delay2','Wire3High','TooEarly'},... %,'Wire3High','TooEarly'
        'OutputActions',{});
    sma = AddState(sma, 'Name','Delay2',...
        'Timer',0,...
        'StateChangeConditions',{'Tup','WaitToPull'},...
        'OutputActions',{'GlobalTimerTrig',1}); % Go signal timing is handled externally; this starts the response window.
    % Column 8 in TrialsMatrix selects the motor/force branch for this trial.
    if S.TrialsMatrix(S.TrialSequence(currentTrial),8)==1
        sma=AddState(sma,'Name', 'WaitToPull',...
            'Timer', 0,...
            'StateChangeConditions',{'GlobalTimer1_End','PreFailAlarm','Wire2High','MotorInterrupt'},...
            'OutputActions', {'WireState',6,'Serial1Code',S.force.thr(currentTrial)});

        sma=AddState(sma,'Name', 'MotorInterrupt',...
            'Timer', 0,...
            'StateChangeConditions',{'GlobalTimer1_End','PreFailAlarm','Wire1High','GracePeriod'},...
            'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});        
        % Require that the joystick remain in the target zone for the configured
        % hold time before a reach is counted as successful.
        sma=AddState(sma,'Name', 'GracePeriod',...
            'Timer', S.GUI.HoldingTime,...
            'StateChangeConditions',{'GlobalTimer1_End', 'PreFailAlarm', 'Tup', 'PreRewardCue','Wire1Low', 'InitMeasurment'},...
            'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});

        sma=AddState(sma,'Name', 'InitMeasurment',...
            'Timer', 0,...
            'StateChangeConditions',{'Tup', 'MotorInterrupt','GlobalTimer1_End', 'PreFailAlarm'},...
            'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});
    else
        sma=AddState(sma,'Name', 'WaitToPull',...
            'Timer', 5,...
            'StateChangeConditions',{'Tup','PreFailAlarm','Wire1High','PreRewardCue'},...
            'OutputActions', {'WireState',10});
    end
   
    sma=AddState(sma,'Name', 'TooEarly',...
        'Timer',0.5,...
        'StateChangeConditions', {'Tup','Break'},...
        'OutputActions', {'SoftCode',1}); % Premature reach: play warning buzzer and skip to the between-trial interval.
    
    sma=AddState(sma,'Name', 'PreRewardCue',...
        'Timer',0.1,...
        'StateChangeConditions', {'Tup','RewardRest'},...
        'OutputActions', {'SoftCode',255});
    sma=AddState(sma,'Name', 'RewardRest',...
        'Timer',0.9,...
         'StateChangeConditions', {'Tup','WaterDeliver'},...
        'OutputActions', {});   %
    sma = AddState(sma, 'Name','WaterDeliver',...
        'Timer',0.04,...
        'StateChangeConditions',{'Tup','Break'},...
        'OutputActions',{'ValveState', S.TrialsMatrix(S.TrialSequence(currentTrial),5)});
    
    %%
    sma=AddState(sma,'Name', 'PreFailAlarm',...
        'Timer',0.1,...
        'StateChangeConditions', {'Tup','FailAlarmCue'},...
        'OutputActions', {'SoftCode',255});  
    sma=AddState(sma,'Name', 'FailAlarmCue',...
        'Timer',0.9,...
        'StateChangeConditions', {'Tup','FailAlarm'},...
        'OutputActions', {});   
    sma=AddState(sma,'Name', 'FailAlarm',...
        'Timer',1,... 
        'StateChangeConditions', {'Tup','Break'},...
        'OutputActions', {'ValveState', S.TrialsMatrix(S.TrialSequence(currentTrial),7),'SoftCode',S.TrialsMatrix(S.TrialSequence(currentTrial),10)});
    %%
    sma=AddState(sma,'Name', 'Break',...
        'Timer',2,...
        'StateChangeConditions', {'Tup','ITI1'},...
        'OutputActions', {'SoftCode',255,'WireState',1}); % Fixed portion of the between-trial interval; also resets output state.
    sma = AddState(sma, 'Name','ITI1',...
        'Timer',iti(currentTrial)-0.5,...
        'StateChangeConditions',{'Tup','ITI2'},...
        'OutputActions',{});
    sma = AddState(sma, 'Name','ITI2',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','exit'},...
        'OutputActions',{});
    SendStateMatrix(sma);
 
%% Start NI-DAQ in the background immediately before the state machine runs

     Nidaq_photometry_Joystick_v3('WaitToStart');
     RawEvents = RunStateMatrix;
    
%% Stop NI-DAQ and cache raw data in BpodSystem.Data
    Nidaq_photometry_Joystick_v3('Stop');
    [Photo470Data,Photo565Data,Traj_X,Traj_Y]=Nidaq_photometry_Joystick_v3('Save');
       if S.GUI.Photometry
        BpodSystem.Data.Nidaq470Data{currentTrial}=Photo470Data; %col1: fluroscence col2: excitation light
        if S.GUI.DbleChanels
        BpodSystem.Data.Nidaq565Data{currentTrial}=Photo565Data;  %col1: fluroscence col2: excitation light
        end
       end 
        BpodSystem.Data.NidaqX{currentTrial}=Traj_X;
        BpodSystem.Data.NidaqY{currentTrial}=Traj_Y;

%% Save trial events, trial metadata, and force settings
if ~isempty(fieldnames(RawEvents))                                          % If trial data was returned
    BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents);            % Computes trial events from raw data
    BpodSystem.Data.TrialSettings = S;
    S.TrialType= Online_TrialType_8class(currentTrial,S.TrialSequence(currentTrial)); % Adds the trial type of the current trial to data
    BpodSystem.Data.TrialSeq(currentTrial) =S.TrialSequence(currentTrial);            % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
    BpodSystem.Data.ResultType(currentTrial) =S.TrialType;   
    BpodSystem.Data.force.thr(currentTrial)=S.force.thr(currentTrial);
    SaveBpodSessionData;                                                    % Saves the field BpodSystem.Data to the current data file
end

%% Update online figures using the saved trial data
[currentOutcome, currentLickEvents]=Online_LickEvents_8class(currentTrial,S.Names.StateToEnd{1,1});
FigLick=Online_LickinJoystick('update',[],[],[],[],[],FigLick,currentTrial,currentOutcome,S.TrialType,currentLickEvents);

%% photometry nidaq 470
[currentNidaq470, nidaqRaw]=Online_NidaqDemod_Joystick(Photo470Data(:,1),nidaq.LED470,S.GUI.LED470_1Freq,S.GUI.LED470_1Amp,S.Names.StateToEnd{1,1},currentTrial,S.GUI.LED470_1Mod);
FigNidaq470=Online_NidaqPlot_8class_470('update',[],[],FigNidaq470,currentNidaq470,nidaqRaw,S.TrialType);
if S.GUI.LED470_2Mod
 [currentNidaq565, nidaqRaw]=Online_NidaqDemod_Joystick(Photo565Data(:,1),nidaq.LED565,S.GUI.LED470_2Freq,S.GUI.LED470_2Amp,S.Names.StateToEnd{1,1},currentTrial,S.GUI.LED470_2Mod);
 FigNidaq565=Online_NidaqPlot_8class_565('update',[],[],FigNidaq565,currentNidaq565,nidaqRaw,S.TrialType);
elseif S.GUI.LED565_Mod
[currentNidaq565, nidaqRaw]=Online_NidaqDemod_Joystick(Photo565Data(:,1),nidaq.LED565,S.GUI.LED565_Freq,S.GUI.LED565_Amp,S.Names.StateToEnd{1,1},currentTrial,S.GUI.LED565_Mod);
 FigNidaq565=Online_NidaqPlot_8class_565('update',[],[],FigNidaq565,currentNidaq565,nidaqRaw,S.TrialType);  
end

 
%% trajectory nidaq
currentNidaq_traj=Online_NidaqPosition(Traj_X,Traj_Y,currentTrial,S.Names.StateToZero{1,1});
FigTraj=Online_Position_8class('update',[],[],FigTraj,currentNidaq_traj,S.TrialType);
 

HandlePauseCondition; % Checks to see if the protocol is paused. IfBpodPath('Ada') so, waits until user resumes.

if BpodSystem.BeingUsed == 0
    return
end
end
end
