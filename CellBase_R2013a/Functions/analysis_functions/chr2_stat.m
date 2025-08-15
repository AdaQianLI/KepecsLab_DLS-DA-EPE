function [E,L,J]= chr2_stat(cellid,Eff,Lat,Jitter)
prs = inputParser;
addRequired(prs,'cellid',@iscellid)

parse(prs,cellid)


E=Eff;
L=Lat;
J=Jitter;
end