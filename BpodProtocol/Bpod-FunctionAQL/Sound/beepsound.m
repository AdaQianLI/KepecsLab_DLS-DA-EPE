function y = beepsound(fs,f,Dur)
% fs = 8192;
% Dur = 0.1;
spaceduration = 0.04;
% f = 800;
nbeeps = 3;
t = linspace(0,Dur,round(Dur*fs));
y = 0.8*sin(2*pi*f*t); % tone
ys = zeros(1,round(spaceduration*fs)); % space
y = [repmat([y ys],[1 nbeeps-1]) y]; % the whole signal
end