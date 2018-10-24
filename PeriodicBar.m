clear all
close all
clc
Rho=2700;
A1=0.02*0.02*pi;
A2=2*A1;
Ee=71e9;
Le=1;

m1=Rho*A1*Le/6; k1=Ee*A1/Le;
m2=Rho*A2*Le/6; k2=Ee*A2/Le;
mc=[2*m1,m1,0;m1,2*m1+2*m2,m2;0,m2,2*m2];
kc=[k1,-k1,0;-k1,k1+k2,-k2;0,-k2,k2];

NC=10; %Number of cells

for ii=1:1001
    freq(ii)=(ii-1)*20;
    KD3=kc-freq(ii)*freq(ii)*mc;
    X21=-KD3(2,1)/KD3(2,2);
    X23=-KD3(2,3)/KD3(2,2);
    KD=[KD3(1,1)+KD3(1,2)*X21, KD3(1,3)+KD3(1,2)*X23; ...
        KD3(3,1)+KD3(3,2)*X21, KD3(3,3)+KD3(3,2)*X23];
    TT=[-KD(1,1)/KD(1,2)                1/KD(1,2)
        KD(2,2)*KD(1,1)/KD(1,2)-KD(2,1) -KD(2,2)/KD(1,2)];
    Lamda(:,ii)=(eig(TT));
    Mew(ii)=acosh(0.5*(Lamda(1,ii)+Lamda(2,ii)));
    %CellAssembly
    KStructure=zeros(NC+1,NC+1);
    BC=zeros(NC+1,1);BC(1)=1;
    CC=zeros(1,NC+1);CC(NC+1)=1;
    for jj=1:NC
      KStructure(jj:jj+1,jj:jj+1)=KStructure(jj:jj+1,jj:jj+1)+KD;
    end
    X1(ii)=20*log([0,1]*inv(KD)*[1;0]);
    XC(ii)=20*log(CC*inv(KStructure)*BC);
end
figure(1); plot(freq,Lamda(1,:),freq,Lamda(2,:)); grid
xlabel('Frequency (rad/sec)')
ylabel('Eigenvalues')
figure(2); plot(freq,real(Mew),freq,abs(imag(Mew))); grid   
xlabel('Frequency (rad/sec)')
ylabel('Propagation Factor')
figure(3); plot(freq,X1); grid   
xlabel('Frequency (rad/sec)')
ylabel('Response Amplitude')
%axis([0,2,-200,200])
figure(4); plot(freq,XC); grid   
xlabel('Frequency (rad/sec)')
ylabel('Response Amplitude')
axis([0,20000,-500,-200])
