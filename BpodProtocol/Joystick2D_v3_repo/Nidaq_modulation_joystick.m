function Modulated_LED=Nidaq_modulation_joystick(amp,freq,mod)
% NIDAQ_MODULATION_JOYSTICK
% Generate the analog output waveform used to modulate an LED during
% photometry acquisition.
%Generates a sin wave for LED amplitude modulation.
global nidaq S

if mod==0 || freq==0 
    Modulated_LED=(amp/2)*ones(nidaq.duration*nidaq.sample_rate,1); % Constant output when modulation is disabled.
else
DeltaT=1/nidaq.sample_rate;
Time=0:DeltaT:(nidaq.duration-DeltaT);
Modulated_LED=amp*(sin(2*pi*freq*Time)+1)/2; % Shifted sine wave in the 0-to-amp range.
Modulated_LED=Modulated_LED';
end
end