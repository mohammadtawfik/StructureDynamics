%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program presents a very simple problem
% of 3-node bar loading
%Written by: Mohammad Tawfik
%Video explaining the code: NONE
%Text about Finite Element Analysis:
% https://www.researchgate.net/publication/321850256_Finite_Element_Analysis_Book_Draft
%Book DOI: 10.13140/RG.2.2.32391.70560
%
%For the Finite Element Course and other courses
% visit http://AcademyOfKnowledge.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing the memory and display
clear all
clc
close all
%Problem Data
NE=2; %number of elements
Length=2.0; %bar length
Width=0.02; %bar width
Thickness=0.02; %bar thickness
Modulus=71e9; %Aluminum modulus of elasticity
Density=2700; %Aluminum density
%Cross-section area
Area=Width*Thickness*pi;
Le=Length/NE; %Element Length
%Element stiffness matrix
Ke=Modulus*Area*[7 -8 1; -8 16 -8;  1 -8 7]/Le/3;
Me=Density*Area*[4 2 -1;  2 16  2; -1  2 4]*Le/30;
%Global stiffness and mass matrix assembly
%Initializing an empty matrix
KGlobal=zeros(2*NE+1,2*NE+1);
MGlobal=zeros(2*NE+1,2*NE+1);
%Assembling the global matrix
for ii=1:2:NE
    IndexE=2*ii-1;
    KGlobal(IndexE:IndexE+2,IndexE:IndexE+2)= ...
                  KGlobal(IndexE:IndexE+2,IndexE:IndexE+2)+Ke;
    MGlobal(IndexE:IndexE+2,IndexE:IndexE+2)= ...
                  MGlobal(IndexE:IndexE+2,IndexE:IndexE+2)+Me;
    KGlobal(IndexE+2:IndexE+4,IndexE+2:IndexE+4)= ...
                  KGlobal(IndexE+2:IndexE+4,IndexE+2:IndexE+4)+2*Ke;
    MGlobal(IndexE+2:IndexE+4,IndexE+2:IndexE+4)= ...
                  MGlobal(IndexE+2:IndexE+4,IndexE+2:IndexE+4)+2*Me;
end
NC=10; %Number of cells

for ii=1:1001
    freq(ii)=(ii-1)*20;
    KD=KGlobal-freq(ii)*freq(ii)*MGlobal;
    KD11=KD(1,1);  KD33=KD(5,5);
    KD13=KD(1,5);  KD31=KD(5,1);
    KD21=KD(2:4,1);KD23=KD(2:4,5);
    KD12=KD(1,2:4);KD32=KD(5,2:4);
    KD22=KD(2:4,2:4);
    X21=-inv(KD22)*KD21;
    X23=-inv(KD22)*KD23;
    KR=[KD11+KD12*X21, KD13+KD12*X23; ...
        KD31+KD32*X21, KD33+KD32*X23];
    TT=[-KR(1,1)/KR(1,2)                1/KR(1,2)
        KR(2,2)*KR(1,1)/KR(1,2)-KR(2,1) -KR(2,2)/KR(1,2)];
    Lamda(:,ii)=(eig(TT));
    Mew(ii)=acosh(0.5*(Lamda(1,ii)+Lamda(2,ii)));
    %CellAssembly
    KStructure=zeros(NC+1,NC+1);
    BC=zeros(NC+1,1);BC(1)=1;
    CC=zeros(1,NC+1);CC(NC+1)=1;
    for jj=1:NC
      KStructure(jj:jj+1,jj:jj+1)=KStructure(jj:jj+1,jj:jj+1)+KR;
    end
    X1(ii)=20*log([0,1]*inv(KR)*[1;0]);
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
