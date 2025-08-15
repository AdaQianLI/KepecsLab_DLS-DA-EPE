function Bpod_Joystick2D_v3_parameter(Param)
%
%
%

global S
    S.Names.Phase={'MotorOFF_Ll_FreeWater','MotorOFF_L2_FixPosition','MotorON_L3_unexpectedF','MotorON_L4_unexpectedF0','MotorOn_L5_3Force','MotorOn_L6_Opto','MotorOn_RwdOmission','MotorOn_RLForce'};
    S.Names.StateToZero={'Delay','blank'};   % {'AirPuff','blank'}
    S.Names.StateToEnd={'Break','blank'};   % {'AirPuff','blank'}    
    S.Names.Rig=Param.rig;
    

%% General Parameters    
    S.GUI.Phase = 1;
    S.GUIMeta.Phase.Style='popupmenu';
    S.GUIMeta.Phase.String=S.Names.Phase;
    S.GUI.MaxTrials=500;
 	S.GUI.Photometry=1;
    S.GUIMeta.Photometry.Style='checkbox';
    S.GUIMeta.Photometry.String='Auto';
    S.GUI.LED565_Mod=0;
    S.GUIMeta.LED565_Mod.Style='checkbox';
    S.GUIMeta.LED565_Mod.String='Auto';
    S.GUI.LED470_2Mod=0;
    S.GUIMeta.LED470_2Mod.Style='checkbox';
    S.GUIMeta.LED470_2Mod.String='Auto';
    S.GUI.LED470_1Mod=1;
    S.GUIMeta.LED470_1Mod.Style='checkbox';
    S.GUIMeta.LED470_1Mod.String='Auto';
    S.GUI.DbleChanels=0;
    S.GUIMeta.DbleChanels.Style='checkbox';
    S.GUIMeta.DbleChanels.String='Auto';
    S.GUI.OptoStim=0;
    S.GUIMeta.OptoStim.Style='checkbox';
    S.GUIMeta.OptoStim.String='Auto';   
    S.GUIPanels.General={'Phase','MaxTrials','Photometry','LED565_Mod','LED470_2Mod','LED470_1Mod','DbleChanels','OptoStim'};
 
%     S.GUI.DelayTime=0.5; % increment in 0.2, till 1.2s
%     S.GUI.ITI=5;
%     S.GUIPanels.Timing={'DelayTime','ITI'};
    
    S.GUI.L0Force=0;
    S.GUI.L1Force=40;
    S.GUI.L2Force=70;
    S.GUI.L3Force=70;
    S.GUI.HoldingTime=0.15;
    S.GUIPanels.Force={'L0Force','L1Force','L2Force','L3Force','HoldingTime'};
    
    S.GUITabs.General={'General','Force'};
    
%% Nidaq and Photometry
	S.GUI.NidaqDuration=45;
    S.GUI.NidaqSamplingRate=6100;
    S.GUI.LED470_1Wavelength=470;
    S.GUI.LED470_1Amp=Param.LED470_1Amp;
    S.GUI.LED470_1Freq=211;
    
    S.GUI.LED470_2Wavelength=470;
    S.GUI.LED470_2Amp=Param.LED470_2Amp;
    S.GUI.LED470_2Freq=331; 
    
    S.GUI.LED565_Wavelength=565;
    S.GUI.LED565_Amp=Param.LED565Amp;
    S.GUI.LED565_Freq=331; 
   
    
    S.GUIPanels.PhotoParam={'NidaqDuration','NidaqSamplingRate',...
                            'LED470_1Wavelength','LED470_1Amp','LED470_1Freq',...
                            'LED470_2Wavelength','LED470_2Amp','LED470_2Freq',...
                            'LED565_Wavelength','LED565_Amp','LED565_Freq'};                        

    S.GUITabs.Photometry={'PhotoParam'};
    
    S.GUI.DecimateFactor_photon=20;
    S.GUI.DecimateFactor_traj=500;
	S.GUI.BaselineBegin=1.5;
    S.GUI.BaselineEnd=2.5;
    S.GUI.NidaqMin=-0.5;
    S.GUI.NidaqMax=1;
    S.GUIPanels.PlotNidaq={'DecimateFactor_photon','DecimateFactor_traj','NidaqMin','NidaqMax','BaselineBegin','BaselineEnd'};
   
    S.GUITabs.PhotoPlot={'PlotNidaq'};
    
    %%
    S.GUI.SoundRamping=0.1;         %sec
    S.GUI.MeanSoundFrequency = 30000;   %Hz
    S.GUI.WidthOfFrequencies=5;
    S.GUI.NumberOfFrequencies=5;
    S.GUI.SinWaveFre1=6000;
    S.GUI.SinWaveFre2=8000;
    S.GUI.BeepFre=500;
    S.GUI.ActionSoundDur=0.3;
    S.GUI.SoundDur=1;
    S.GUI.SF=192000; % Sound card sampling 
    S.GUI.LowFreq=5000;
    S.GUI.HighFreq=40000;
    S.GUI.DelayDur=18;
    S.GUIPanels.SoundParameter={'SoundRamping','MeanSoundFrequency','WidthOfFrequencies','NumberOfFrequencies','SinWaveFre1','SinWaveFre2','BeepFre','ActionSoundDur','SoundDur','SF','LowFreq','HighFreq','DelayDur'};

    S.GUITabs.Sound={'SoundParameter'};
    
%% Online Plots   
    S.GUI.StateToZero=1;
	S.GUIMeta.StateToZero.Style='popupmenu';
    S.GUIMeta.StateToZero.String=S.Names.StateToZero;
    S.GUI.StateToEnd=1;
	S.GUIMeta.StateToEnd.Style='popupmenu';
    S.GUIMeta.StateToEnd.String=S.Names.StateToEnd;
    S.GUI.TimeMin_lick=-5;
    S.GUI.TimeMax_lick=3;
    S.GUI.TimeMin_photometry=-5;
    S.GUI.TimeMax_photometry=3;
    S.GUIPanels.PlotParameters={'StateToZero','StateToEnd','TimeMin_lick','TimeMax_lick','TimeMin_photometry','TimeMin_photometry'};
       
    S.GUITabs.OnlinePlot={'PlotParameters'};
end
