function y = beepsound(fs,f,Dur)
% BEEPSOUND
% Generate a short multi-beep warning tone. In this protocol it is used for
% early reaches before the go signal (SoftCode 4).
% fs = 8192;
% Dur = 0.1;
spaceduration = 0.04; % Gap between beeps (s).
% f = 800;
nbeeps = 3; % Three brief warning beeps.
t = linspace(0,Dur,round(Dur*fs));
y = 0.8*sin(2*pi*f*t); % tone
ys = zeros(1,round(spaceduration*fs)); % space
y = [repmat([y ys],[1 nbeeps-1]) y]; % Full warning waveform used by SoftCode 4.
end