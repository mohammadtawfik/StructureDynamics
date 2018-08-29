clear all
close all

Es=1;
Ep=1;
Cd=1;

for ii=1:101
  Omega=(ii-1)*.05;
  Freq(ii)=Omega;
  EComplex=(Es*Ep+Es*Cd*i*Omega)/(Es+Ep+(i*Omega*Cd))
  Estorage(ii)=real(EComplex);
  Eloss=imag(EComplex);
  Eta(ii)=Eloss/Estorage(ii);
end

plot(Freq,Eta,Freq,Estorage)
axis([0,5,0,1])
grid