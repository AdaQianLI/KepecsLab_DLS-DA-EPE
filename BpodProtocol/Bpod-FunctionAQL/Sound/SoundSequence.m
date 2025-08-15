function Seq=SoundSequence(Reps,NbFreq)
    List=1:1:NbFreq;
    List=List(randperm(length(List)));
    Seq=zeros(NbFreq*Reps,1);
    for i=1:NbFreq
        Seq((Reps*(i-1))+1:(Reps*(i)))=List(i);
    end
end