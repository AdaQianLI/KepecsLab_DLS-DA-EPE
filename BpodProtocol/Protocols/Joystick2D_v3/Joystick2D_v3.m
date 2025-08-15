function Joystick2D_v3
global BpodSystem nidaq S

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
ParamPC=BpodParam_PCdep_Joystick_v3();
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    Bpod_Joystick2D_v3_parameter(ParamPC);
end

% Initialize parameter GUI plugin and Pause
BpodParameterGUI('init', S);
BpodSystem.Pause=1;
HandlePauseCondition;
 S = BpodParameterGUI('sync', S);
 
%% Define Sound and send to sound serve
FreeWaterCue=GenerateBeatWave(S.GUI.SF,S.GUI.SinWaveFre2,10);
TimeOutCue = WhiteNoise(S.GUI.SF,S.GUI.SoundDur);  %noise
SuccCue=GenerateBeatWave1(S.GUI.SF,S.GUI.SinWaveFre1,S.GUI.SoundDur);
TooEealyCue=beepsound(S.GUI.SF,S.GUI.BeepFre,0.1);
% Program sound server
PsychToolboxSoundServer('init');
PsychToolboxSoundServer('Load', 1, FreeWaterCue);
PsychToolboxSoundServer('Load', 2, TimeOutCue);  %noise 
PsychToolboxSoundServer('Load', 3, SuccCue); 
PsychToolboxSoundServer('Load', 4, TooEealyCue); 
BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler_PlaySound';


 [S.TrialsMatrix, S.TrialsNames, S.TrialsResNames]=Joystick2D_v3_Phase(S.Names.Phase{S.GUI.Phase});  % choose Air_Stim_Pairing
 tempSequence=WeightedRandomTrials(S.TrialsMatrix(:,2)', S.GUI.MaxTrials);

%  idx=find(tempSequence==1);
%  idx_=idx+1;
%  temp=find(tempSequence(idx_)==1);
%  tempSequence(idx_(temp)-1)=3;
 S.TrialSequence=[tempSequence];%[ones(1,100) ones(1,100)*2 ones(1,100)*3 ones(1,100)*4];  %ones(1,30) 
 
for ii=1:length(S.TrialSequence)
    switch S.TrialSequence(ii)
        case 1
            S.force.thr(ii)=S.GUI.L0Force;%round((S.GUI.L2Force-S.GUI.L1Force)*rand(1)+S.GUI.L1Force);
        case 2
            S.force.thr(ii)=S.GUI.L1Force;
        case 3
            S.force.thr(ii)=S.GUI.L2Force;%round((S.GUI.L3Force-S.GUI.L2Force)*rand(1)+S.GUI.L2Force);    
        case 4
            S.force.thr(ii)=S.GUI.L3Force;%opto;    
    end
end
if S.GUI.OptoStim == 1
    PulsePal('COM5');
    load('C:\Users\Kepecs\Documents\Data\Ada\Bpod\Protocols\Joystick2D_v3_noRwdcueNew\LightTrain_500ms.mat');
    ProgramPulsePal(ParameterMatrix);
end
%% ITI
n=length(S.TrialSequence);
lambda=1/3;
a=1;
b=6;
fa=1-exp(-lambda*a);
fb=1-exp(-lambda*b);
u=fa+(fb-fa)*rand(n,1);
iti=-log(1-u)/lambda;
iti=iti';

%% NIDAQ Initialization

Nidaq_photometry_Joystick_v3('ini',ParamPC);
FigLick=Online_LickinJoystick('ini',S.TrialSequence,S.TrialsMatrix,S.Names.Phase{S.GUI.Phase},S.TrialsNames,S.TrialsResNames);
FigNidaq470=Online_NidaqPlot_8class_470('ini',S.Names.Phase{S.GUI.Phase},S.TrialsResNames); %'Air_Stim_Pairing'
if S.GUI.DbleChanels == 1
FigNidaq565=Online_NidaqPlot_8class_565('ini',S.Names.Phase{S.GUI.Phase},S.TrialsResNames); %'Air_Stim_Pairing'
end    
FigTraj=Online_Position_8class('ini',S.Names.Phase{S.GUI.Phase},S.TrialsResNames);


%% Main trial loop
for currentTrial = 1:length(S.TrialSequence)
S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin 
      
%% Assemble State matrix
 	sma = NewStateMatrix();
    sma = SetGlobalTimer(sma,1,5.5);
    sma = SetGlobalTimer(sma,3,S.GUI.HoldingTime);
    sma = AddState(sma, 'Name','Dummy1',...
        'Timer',4,...  %4
        'StateChangeConditions',{'Tup','CueDeliver'},...
        'OutputActions',{}); 
    sma = AddState(sma, 'Name','CueDeliver',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','Delay','Wire3High','TooEarly'},...
        'OutputActions',{'BNCState',S.TrialsMatrix(S.TrialSequence(currentTrial),9)});   
    sma = AddState(sma, 'Name','Delay',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','Delay1','Wire3High','TooEarly'},...
        'OutputActions',{});
    sma = AddState(sma, 'Name','Delay1',...
        'Timer',0.5,...
        'StateChangeConditions',{'Tup','Delay2','Wire3High','TooEarly'},...
        'OutputActions',{'BNCState',S.TrialsMatrix(S.TrialSequence(currentTrial),3)});
    sma = AddState(sma, 'Name','Delay2',...
        'Timer',0,...
        'StateChangeConditions',{'Tup','WaitToPull'},...
        'OutputActions',{'GlobalTimerTrig',1});
    if S.TrialsMatrix(S.TrialSequence(currentTrial),8)==1
        sma=AddState(sma,'Name', 'WaitToPull',...
            'Timer', 0,...
            'StateChangeConditions',{'GlobalTimer1_End','PreFailAlarm','Wire2High','MotorInterrupt'},...
            'OutputActions', {'WireState',6,'Serial1Code',S.force.thr(currentTrial)});
%         sma=AddState(sma,'Name', 'ForceWithhold',...
%             'Timer',0.5+rand*0.2,...%+rand*0.3,...%0.4+rand*0.3,..
%             'StateChangeConditions',{'GlobalTimer1_End','PreFailAlarm','Tup','MotorInterrupt'},...
%             'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});  
        sma=AddState(sma,'Name', 'MotorInterrupt',...
            'Timer', 0,...
            'StateChangeConditions',{'GlobalTimer1_End','PreFailAlarm','Wire1High','InitMeasurment'},...
            'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Measurement of time the subject kept the joystick in the right place %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sma=AddState(sma,'Name', 'InitMeasurment',...
            'Timer', 0,...
            'StateChangeConditions',{'GlobalTimer1_End', 'PreFailAlarm', 'Tup', 'MotorInterrupt2'},...
            'OutputActions', {'GlobalTimerTrig', 3, 'WireState',10,'Serial1Code',S.force.thr(currentTrial)});
        sma=AddState(sma,'Name', 'MotorInterrupt2',...
            'Timer', 0,...
            'StateChangeConditions',{'GlobalTimer3_End', 'PreRewardCue', 'Wire1Low', 'GracePeriod', 'GlobalTimer1_End', 'PreFailAlarm'},...
            'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});
        sma=AddState(sma,'Name', 'GracePeriod',...
            'Timer', 0.005,...
            'StateChangeConditions',{'Wire1High', 'MotorInterrupt2', 'Tup', 'MotorInterrupt', 'GlobalTimer3_End','MotorInterrupt', 'GlobalTimer1_End', 'PreFailAlarm'},...
            'OutputActions', {'WireState',10,'Serial1Code',S.force.thr(currentTrial)});
    else
        sma=AddState(sma,'Name', 'WaitToPull',...
            'Timer', 3,...
            'StateChangeConditions',{'Tup','PreRewardCue','Wire1High','PreRewardCue'},...
            'OutputActions', {'WireState',10,'SoftCode',1});
    end
   
    sma=AddState(sma,'Name', 'TooEarly',...
        'Timer',0.5,...
        'StateChangeConditions', {'Tup','Break'},...
        'OutputActions', {'SoftCode',4});
    
    sma=AddState(sma,'Name', 'PreRewardCue',...
        'Timer',0.1,...
        'StateChangeConditions', {'Tup','RewardRest'},...
        'OutputActions', {'SoftCode',255});
    sma=AddState(sma,'Name', 'RewardRest',...
        'Timer',0.9,...    %1
        'StateChangeConditions', {'Tup','WaterDeliver'},...
        'OutputActions', {});   %
    sma = AddState(sma, 'Name','WaterDeliver',...
        'Timer',0.1,...%S.TrialsMatrix(S.TrialSequence(currentTrial),11),...
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
        'Timer',0.5,...
        'StateChangeConditions', {'Tup','Break'},...
        'OutputActions', {'ValveState', S.TrialsMatrix(S.TrialSequence(currentTrial),7),'SoftCode',2});
    %%
    sma=AddState(sma,'Name', 'Break',...
        'Timer',2,...
        'StateChangeConditions', {'Tup','ITI'},...
        'OutputActions', {'SoftCode',255});
    sma = AddState(sma, 'Name','ITI',...
        'Timer',iti(currentTrial),...
        'StateChangeConditions',{'Tup','exit'},...
        'OutputActions',{'WireState',1});
    SendStateMatrix(sma);
 
%% NIDAQ Get nidaq ready to start

     Nidaq_photometry_Joystick_v3('WaitToStart');
     RawEvents = RunStateMatrix;
    
%% NIDAQ Stop acquisition and save data in bpod structure
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
%     BpodSystem.Data.NidaqData{currentTrial} = nidaq.ai_data;
%% Save
if ~isempty(fieldnames(RawEvents))                                          % If trial data was returned
    BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents);            % Computes trial events from raw data
    BpodSystem.Data.TrialSettings = S;
    S.TrialType= Online_TrialType_8class(currentTrial,S.TrialSequence(currentTrial)); % Adds the trial type of the current trial to data
    BpodSystem.Data.TrialSeq(currentTrial) =S.TrialSequence(currentTrial);                       % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
    BpodSystem.Data.ResultType(currentTrial) =S.TrialType;   
    BpodSystem.Data.force.thr(currentTrial)=S.force.thr(currentTrial);
    SaveBpodSessionData;                                                    % Saves the field BpodSystem.Data to the current data file
end

%% PLOT - extract events from BpodSystem.data and update figures
%  S.TrialType= Online_TrialType_8class(currentTrial,S.TrialSequence(currentTrial)); % Adds the trial type of the current trial to data
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
% currentNidaq_470=Online_NidaqEvents(nidaq.ai_data(:,1),currentTrial,S.Names.StateToEnd{1,1});

% FigNidaq=Online_NidaqPlot_8class('update',[],FigNidaq,currentNidaq_470,S.TrialType);

%% photometry nidaq 565
% if S.GUI.DbleChanels
%     [currentNidaq565,nidaqRawb]=Online_NidaqDemod_Joystick(Photo565Data(:,1),nidaq.LED565,S.GUI.LED565_Freq,S.GUI.LED565_Amp,S.Names.StateToEnd{1,1},currentTrial);
%     FigNidaq565=Online_NidaqPlot_8class_565('update',[],FigNidaq565,currentNidaq565,nidaqRawb,S.TrialType);
% end
 
%% trajectory nidaq
currentNidaq_traj=Online_NidaqPosition(Traj_X,Traj_Y,currentTrial,S.Names.StateToZero{1,1});
FigTraj=Online_Position_8class('update',[],[],FigTraj,currentNidaq_traj,S.TrialType);
 

HandlePauseCondition; % Checks to see if the protocol is paused. IfBpodPath('Ada') so, waits until user resumes.

if BpodSystem.BeingUsed == 0
    return
end
end
end
