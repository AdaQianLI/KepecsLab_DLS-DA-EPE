function [events,epochs] = defineEventsEpochs_Joystick2D
%DEFINEEVENTSEPOCHS_DEFAULT   Define events and epochs for spike extraction.
%   [EVENTS,EPOCHS] = DEFINEEVENTSEPOCHS_DEFAULT defines events and epochs
%   for spike extraction. 
%
%   EVENTS is a Nx4 cell array with columns corresponding to EventLabel,
%   EventTrigger1, EventTrigger2, Window. EventLabel is the name for
%   referencing the event. EventTrigger1 and EventTrigger2 are names of
%   TrialEvent variables (e.g. 'LeftPortIn'). For fixed windows, the two
%   events are the same; for variable windows, they correspond to the start
%   and end events. Window specifies time offsets relative to the events;
%   e.g. events(1,:) = {'OdorValveOn','OdorValveOn','OdorValveOn',[-3 3]};
%
%   EPOCH is a Nx4 cell array with columns corresponding to  EpochLabel, 
%   ReferenceEvent, Window, RealWindow. EventLabel is the name for 
%   referencing the epoch. ReferenceEvent should match an EventLabel in 
%   EVENTS (used for calculating the epoch rates). RealWindow is currently
%   not implemented (allocated for later versions).
%
%   See also MAKETRIALEVENTS2_GONOGO and DEFINEEVENTSEPOCHS_PULSEON.

%   Edit log: BH 7/6/12

% Define events and epochs
%              EventLabel       EventTrigger1      EventTrigger2      Window
i = 1;
events(i,:) = {'AccemaxTime',       'AccemaxTime',         'AccemaxTime',         [-3 3]};    i = i + 1;
events(i,:) = {'AcceXmaxTime',      'AcceXmaxTime',        'AcceXmaxTime',        [-3 3]};    i = i + 1;
events(i,:) = {'Break',             'Break',               'Break',               [-3 3]};    i = i + 1;
events(i,:) = {'CueDeliver',        'CueDeliver',          'CueDeliver',          [-3 3]};    i = i + 1;
   
events(i,:) = {'Delay',             'Delay',               'Delay',                [-3 3]};    i = i + 1;
events(i,:) = {'Delay2',            'Delay2',              'Delay2',               [-3 3]};    i = i + 1;
events(i,:) = {'Dummy1',            'Dummy1',              'Dummy1',               [-3 3]};    i = i + 1;
events(i,:) = {'FailAlarm',         'FailAlarm',           'FailAlarm',            [-3 3]};    i = i + 1;
events(i,:) = {'FailAlarmCue',       'FailAlarmCue',       'FailAlarmCue',         [-3 3]};    i = i + 1;
events(i,:) = {'ITI',                'ITI',                'ITI',                  [-3 3]};    i = i + 1;
events(i,:) = {'MotorInterrupt',     'MotorInterrupt',     'MotorInterrupt',       [-3 3]};    i = i + 1;

events(i,:) = {'OffsetTime',         'OffsetTime',         'OffsetTime',           [-3 3]};    i = i + 1;
events(i,:) = {'OnsetTime',          'OnsetTime',          'OnsetTime',            [-3 3]};   i = i + 1;
events(i,:) = {'BurstPeriod',        'PreFailAlarm',       'PreFailAlarm',         [-3 3]};   i = i + 1;
events(i,:) = {'PreRewardCue',       'PreRewardCue',       'PreRewardCue',         [-3 3]};    i = i + 1;
events(i,:) = {'VelmaxTime',         'VelmaxTime',         'VelmaxTime',           [-3 3]};    i = i + 1;
events(i,:) = {'VelXmaxTime',        'VelXmaxTime',        'VelXmaxTime',          [-3 3]};    i = i + 1;
events(i,:) = {'WaitToPull',         'WaitToPull',         'WaitToPull',           [-3 3]};    i = i + 1;
events(i,:) = {'WaterDeliver',       'WaterDeliver',       'WaterDeliver',         [-3 3]};    i = i + 1;
events(i,:) = {'XmaxTime',           'XmaxTime',           'XmaxTime',             [-3 3]};    i = i + 1;
% Define epochs for rate calculations
%               EpochLabel      ReferenceEvent     Window             RealWindow
 i = 1;
 
epochs(i,:) = {'Dummy1',    'Dummy1',       [0 3],         'StimulusSampling'};    i = i + 1;
epochs(i,:) = {'Delay',     'Delay',        [0 1],         'StimulusSampling'};    i = i + 1;
% epochs(i,:) = {'FixedLightResponse','PulseOn',     [0.0 0.005],        'PulsePeriod'};    i = i + 1;
% 
% epochs(i,:) = {'Baseline',       'PrePulseIPI',    [NaN NaN],          'NaN'};            i = i + 1;
% epochs(i,:) = {'LightResponse',  'PulsePeriod',    [NaN NaN],          'NaN'};            i = i + 1;