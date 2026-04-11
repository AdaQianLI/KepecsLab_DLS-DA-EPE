function Sound = GenerateBeatWave1(SamplingRate, Frequency, Duration)
% GENERATEBEATWAVE1
% Build the pulsed beat-wave used here as the success cue. The waveform
% contains three tone packets separated by silent gaps.
% Duration in seconds
dt = 1/SamplingRate;

t = 0:dt:Duration/5; % One pulse occupies one fifth of the requested overall duration.
FreqDiff = 5;
factor=(Frequency/25000).^2;
BeatWave=(10*sin(2*pi*Frequency*t) + sin(2*pi*(Frequency-FreqDiff)*t))*factor;
iti=zeros(1,Duration/5*SamplingRate); % Silent gap between pulses.
Sound=[BeatWave iti BeatWave iti BeatWave]; % Three-pulse success cue.

