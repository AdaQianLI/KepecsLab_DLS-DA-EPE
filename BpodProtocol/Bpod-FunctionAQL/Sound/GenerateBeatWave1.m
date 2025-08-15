function Sound0 = GenerateBeatWave1(SamplingRate, Frequency, Duration)
% Duration in seconds
dt = 1/SamplingRate;

t = 0:dt:Duration/5;
FreqDiff = 5;
factor=(Frequency/25000).^2;
BeatWave=(10*sin(2*pi*Frequency*t) + sin(2*pi*(Frequency-FreqDiff)*t))*factor;
iti=zeros(1,Duration/5*SamplingRate);
Sound0=[BeatWave iti BeatWave iti BeatWave];
end
