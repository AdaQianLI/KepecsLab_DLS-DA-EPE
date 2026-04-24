function Param=BpodParam_PCdep_Joystick_v3()
% BPODPARAM_PCDEP_JOYSTICK_V3
% Return rig-specific defaults that depend on the acquisition computer.
% This keeps NI-DAQ device names and default LED amplitudes out of the
% main protocol file.

% Edit or extend this switch block when moving the protocol to a new acquisition PC.
switch getenv('computername')
    case 'KEPECSPHOTO-01'
        Param.rig='Photometry1';
        Param.nidaqDev='Dev1';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.3;
        Param.LED565Amp=0.5;
    case 'KEPECSPHOTO-02'
        Param.rig='Photometry2';
        Param.nidaqDev='Dev3';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.3;
        Param.LED565Amp=0.5;
    case 'KEPECSPHOTO-05' %mid rig
        Param.rig='Photometry5';
        Param.nidaqDev='Dev5';
        Param.LED470_1Amp=0.18;
        Param.LED470_2Amp=0;
        Param.LED565Amp=1.2;
	case 'KEPECSPHOTO-06'
        Param.LED565Amp=5;
        Param.rig='Photometry6';
        Param.nidaqDev='Dev2';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.2;
        Param.LED565Amp=0.4;
end
end