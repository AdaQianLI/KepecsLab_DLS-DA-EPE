function Param=BpodParam_PCdep()

switch getenv('computername')
    case 'KEPECSPHOTO-01'
        Param.rig='Photometry1';
        Param.nidaqDev='Dev2';
        Param.LED1Amp=0.3;
        Param.LED2Amp=5;
        Param.LED1bAmp=0;
    case 'KEPECSPHOTO02'
        Param.rig='Photometry2';
        Param.nidaqDev='Dev2';
        Param.LED1Amp=0.3;
        Param.LED2Amp=5;
        Param.LED1bAmp=0;
    case 'KEPECSPHOTO-03'
        Param.rig='Photometry3';
        Param.nidaqDev='Dev2';
        Param.LED1Amp=0.5;
        Param.LED2Amp=5;
        Param.LED1bAmp=0;
	case 'KEPECSPHOTO-06'
        Param.rig='Photometry6';
        Param.nidaqDev='Dev2';
        Param.LED1Amp=0.3;
        Param.LED2Amp=5;
        Param.LED1bAmp=0;
end
end