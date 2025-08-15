function [CueC,CueD]=SoundsOfTheDay(DayNb,weekNB)

%% List of available numbers said by different people
% Week1 - MEHA (zero 1 - 2 3 - 4 5 - 6 7 - 8 9)
switch weekNB
    case 1
ListOfSound=...
    {'FEA_ZA.wav','FEA_2A.wav','FEA_4A.wav','FEA_6A.wav','FEA_8A.wav';...
     'FEA_1A.wav','FEA_3A.wav','FEA_5A.wav','FEA_7A.wav','FEA_9A.wav'};
    case 2
ListOfSound=...
    {'MEH_ZA.wav','MEH_2A.wav','MEH_4A.wav','MEH_6A.wav','MEH_8A.wav';...
     'MEH_1A.wav','MEH_3A.wav','MEH_5A.wav','MEH_7A.wav','MEH_9A.wav'}; 

end

%% Open audio files
cd 'C:\Users\Kepecs\Documents\Data\Quentin\Function QC\NumberSounds'
CueC=audioread(ListOfSound{1,DayNb});
CueC=resample(CueC',192000,8000)*300;
CueD=audioread(ListOfSound{2,DayNb});
CueD=resample(CueD',192000,8000)*300;

end