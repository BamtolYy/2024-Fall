function [jk,Nk]        = calculateSampleBounds(tsk,fs,Taccumulation,Tsample);
jk = round(tsk(mm).*fs)+1;
Nk = Taccumulation/Tsample;
if Nk