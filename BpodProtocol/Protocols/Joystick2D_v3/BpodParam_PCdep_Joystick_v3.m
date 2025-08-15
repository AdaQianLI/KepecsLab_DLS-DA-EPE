function Param=BpodParam_PCdep_Joystick_v3()

switch getenv('computername')
    case 'KEPECSPHOTO-01'
        Param.rig='Photometry1';
        Param.nidaqDev='Dev1';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.3;
        Param.LED565Amp=5;
    case 'KEPECSPHOTO02'
        Param.rig='Photometry2';
        Param.nidaqDev='Dev3';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.3;
        Param.LED565Amp=5;
    case 'KEPECSPHOTO-03'
        Param.rig='Photometry3';
        Param.nidaqDev='Dev1';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.3;
	case 'KEPECSPHOTO-06'
        Param.LED565Amp=5;
        Param.rig='Photometry6';
        Param.nidaqDev='Dev2';
        Param.LED470_1Amp=0.3;
        Param.LED470_2Amp=0.2;
        Param.LED565Amp=0.4;
end
end