function SoftCodeHandler_PlaySound(SoundID)
% SOFTCODEHANDLER_PLAYSOUND
% Play or stop PsychToolbox sounds in response to Bpod soft codes.
% SoundID 255 is reserved here as a stop-all command.
if SoundID ~= 255
    PsychToolboxSoundServer('Play', SoundID);
else
    PsychToolboxSoundServer('StopAll');
end