function trialSeq=WeightedRandomTrials(probTrials,maxTrials)
%trialSeq=WeightedRandomTrials(probTrials,maxTrials).
%
%Generates a randomized and weighted distribution "TrialSeq"
%The values of the random distribution are of the size of "probTrials".
%"proTrials" defines the occurence probabilities of the different values
%"maxTrials" defines the size of the random distribution
%
%function written by Quentin for DelayedReward bpod protocol

if nargin~=2
    disp('*** please enter 2 arguments for the WeightedRandomTrials function ***');
    return
end
if abs(sum(probTrials))-1 > 1e-9
    disp('*** Error in defineRandomizedTrials, typeMatrix proportions do not add up to 1 ***');
    trialSeq = [];
    return
end

rng('shuffle')
tempSeq=rand(1,maxTrials);
trialSeq=arrayfun(@(z)sum(z>=cumsum([0,probTrials])),tempSeq);
end