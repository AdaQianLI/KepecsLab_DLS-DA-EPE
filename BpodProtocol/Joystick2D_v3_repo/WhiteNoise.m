function sound=WhiteNoise(sampRate,duration)
% WHITENOISE
% Generate white noise for the timeout/alarm sound.
%sound=SoundGenerator(sampRate, meanFreq, widthFreq, nbOfFreq, duration,rampTime).
%
%Generates a ramping sound constituted by multiple frequencies.
%The frequencies are defined by "meanFreq", "widthFreq" and "nbOfFreq".
%SamplingRate is the sampling Rate of the sound card.
%
%function written by Quentin for DelayedReward bpod protocol

dt = 1/sampRate;
t = 0:dt:duration;
sound=rand(1,length(t)); % Uniform white-noise segment.
