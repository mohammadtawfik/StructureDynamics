clear all
close all

Es=1;
Cd=1;

for ii=1:101
  Omega=(ii-1)*.05;
  Freq(ii)=Omega;
  Estorage(ii)=Es;
  Eloss=Cd*Omega;
  Eta(ii)=Eloss/Estorage(ii);
end

plot(Freq,Eta,Freq,Estorage)
axis([0,5,0,2])
grid
