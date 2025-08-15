function sma=ScaledRewardStates(sma, firstState, postState, gTimer, nbOfStates, minValveTime, maxValveTime, rewardValveCode, BufTime)
%sma=ScaledRewardStates(sma, firstState, postState, gTimer, nbOfStates,minValveTime, maxValveTime, rewardValveCode, BufTime).
%Add states to a preexisting state matrix to deliver a ramping reward.
%The previous state need to be "firstState" and the following state need to be "postState".
%This function introduce a delay "gTimer" between the firstState - usually the stimulus.
%GlobalTimer #5 computes the delay, GlobaTimer #4 fixes the time between the reward and the post reward State.
%Reward amounts are discretized by a "nbOfStates" from a "minValveTime" to "maxValveTime" amount, and delivered through the "rewardValveCode".
%function written by Quentin for DelayedReward bpod protocol

if nargin~=9
    disp('*** please enter 9 arguments for the ScaledRewardStates function ***');
    return
end

reward_valve_times = linspace(minValveTime,maxValveTime,nbOfStates);
timeStep=gTimer/nbOfStates;

sma = SetGlobalTimer(sma, 5, gTimer);                  % to normalize time between cue and reward
sma = SetGlobalTimer(sma, 4, maxValveTime+BufTime);    % to normalize time between reward and poststate

sma = AddState(sma, 'Name', firstState, ...
    'Timer', 0,...
    'StateChangeConditions', {'Tup','WaitForLick_1'},...
    'OutputActions', {'GlobalTimerTrig',5});

for i = 1:nbOfStates
    % The main path - if the animal licks.
    this_wait_state         = sprintf('WaitForLick_%.0f',i);
    this_global_timer_state = sprintf('GlobalTimer_%.0f',i);
    this_reward_state       = sprintf('RewardSate_%.0f',i);

    % What to do at the end of this state?
    if i < nbOfStates,
            next_wait_state = sprintf('WaitForLick_%.0f',i+1);
    else
        next_wait_state = 'NeverLickedState';
    end

    this_reward_valve_time = reward_valve_times(i);

    % Wait for lick
    sma = AddState(sma, 'Name', this_wait_state, ...
        'Timer',timeStep,...
        'StateChangeConditions', {'Port2In',this_global_timer_state,'Tup',next_wait_state},...
        'OutputActions', {});
    
    % If lick, use this global timer.
    sma = AddState(sma,'Name', this_global_timer_state, ...
        'Timer',20,...
        'StateChangeConditions', {'GlobalTimer5_End', this_reward_state},...
        'OutputActions', {});

    % If licked, use this reward
    sma = AddState(sma,'Name', this_reward_state, ...
        'Timer',this_reward_valve_time,...
        'StateChangeConditions', {'Tup', 'RewardWaitTimer'},...
        'OutputActions', {'ValveState', rewardValveCode,'GlobalTimerTrig',4});   
end
    % If no lick
    sma=AddState(sma,'Name', 'NeverLickedState',...
        'Timer', 0,...
        'StateChangeConditions', {'Tup', 'RewardWaitTimer'},...
        'OutputActions', {'GlobalTimerTrig',4});
    
    sma=AddState(sma,'Name', 'RewardWaitTimer',...
        'Timer', 0,...
        'StateChangeConditions', {'GlobalTimer4_End', postState},...
        'OutputActions', {});
end