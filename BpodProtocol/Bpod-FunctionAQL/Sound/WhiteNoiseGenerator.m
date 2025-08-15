function noise = WhiteNoiseGenerator(sampRate,duration,rampTime)

if nargin==2
    rampTime=0;
end
    
SoundLength=duration*sampRate;
noise=2*rand(SoundLength,1)-1;

ampl=ones(1,duration*sampRate)*20;
ampl(1,1:rampTime*sampRate)=linspace(0,1,sampRate*rampTime);
ampl(1,sampRate*(duration-rampTime)+1:end)=linspace(1,0,sampRate*rampTime);

noise=ampl.*noise';
    
end