clear all
close all
clc

m=1; k=1;
mc=[m,0,0;0,m,0;0,0,m];
kc=[k,-k,0;-k,2*k,-k;0,-k,k];
M=2*m;
K=2*k;
M7=[M,0,0,0,0,0,0; ...
    0,m,0,0,0,0,0; ...
    0,0,M,0,0,0,0; ...
    0,0,0,m,0,0,0; ...
    0,0,0,0,M,0,0; ...
    0,0,0,0,0,m,0; ...
    0,0,0,0,0,0,M];
K7=[k,-k,0,0,0,0,0; ...
    -k,K,-k,0,0,0,0; ...
    0,-k,K,-k,0,0,0; ...
    0,0,-k,K,-k,0,0; ...
    0,0,0,-k,K,-k,0; ...
    0,0,0,0,-k,K,-k; ...
    0,0,0,0,0,-k,k];
for ii=1:1001
    freq(ii)=(ii-1)*0.003;
    KD3=kc-freq(ii)*freq(ii)*mc;
    KD7=K7-freq(ii)*freq(ii)*M7;
    X3(ii)=20*log([0,0,1]*pinv(KD3)*[1;0;0]);
    X7(ii)=20*log([0,0,0,0,0,0,1]*pinv(KD7)*[1;0;0;0;0;0;0]);
    X21=-KD3(2,1)/KD3(2,2);
    X23=-KD3(2,3)/KD3(2,2);
    KD=[KD3(1,1)+KD3(1,2)*X21, KD3(1,3)+KD3(1,2)*X23; ...
        KD3(3,1)+KD3(3,2)*X21, KD3(3,3)+KD3(3,2)*X23];
    TT=[-KD(1,1)/KD(1,2)                1/KD(1,2)
        KD(2,2)*KD(1,1)/KD(1,2)-KD(2,1) -KD(2,2)/KD(1,2)];
    Lamda(:,ii)=(eig(TT));
    Mew(ii)=acosh(0.5*(Lamda(1,ii)+Lamda(2,ii)));
end
figure(1); plot(freq,Lamda(1,:),freq,Lamda(2,:)); grid
xlabel('Frequency (rad/sec)')
ylabel('Eigenvalues')
axis([0,2,-5,5])
figure(2); plot(freq,real(Mew),freq,abs(imag(Mew))); grid   
xlabel('Frequency (rad/sec)')
ylabel('Propagation Factor')
axis([0,2,0,3.5])
figure(3); plot(freq,X3); grid   
xlabel('Frequency (rad/sec)')
ylabel('Response Amplitude')
axis([0,2,-200,200])
figure(4); plot(freq,X7); grid   
xlabel('Frequency (rad/sec)')
ylabel('Response Amplitude')
axis([0,2,-200,200])
