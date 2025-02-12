%% PIRD for N=2 source processes - state space computation

% it -  index of the target process (1,2 or 3)

function out=pird2(A,C,K,Su,Mv,it,Fs,nfft)

narginchk(6,8)
if nargin<8, nfft=512; end
if nargin<7, Fs=1; end %default non compute time domain measures through submodels

is=setdiff(1:3,it); %index of source processes

outtmp=oir_mir(A,C,K,Su,Mv,it,is(1),Fs,nfft);
i1=0.5*outtmp.f12; % spectral MIR between Y and X1
I1=mean(i1); % time-domain MIR between Y and X1
outtmp=oir_mir(A,C,K,Su,Mv,it,is(2),Fs,nfft);
i2=0.5*outtmp.f12; % spectral MIR between Y and X2
I2=mean(i2); % time-domain MIR between Y and X2

outtmp=oir_mir(A,C,K,Su,Mv,it,is,Fs,nfft);
i12=0.5*outtmp.f12; % spectral MIR between Y and [X1,X2]
I12=mean(i12); % time-domain MIR between Y and [X1,X2]

% Spectral OIR
outtmp=oir_deltaO(A,C,K,Su,Mv,1:3,it);
o123=0.5*outtmp.dO12f;
% iir=i1+i2-i12; % verification: spectral interaction information rate

% Redundancy functions
icap_min=min([i1 i2],[],2); %spectral redundancy based on local minimum MI in frequency
icap_oi=o123.*(o123>0); %spectral redundancy based on local OIR in frequency - not used (can lead to negative atmos - to be explored)

Icap_min=mean(icap_min);
Icap_oi=mean(icap_oi);

R_min=Icap_min;
U1_min=I1-R_min; U2_min=I2-R_min;
S_min=I12-R_min-U1_min-U2_min;

R_oi=Icap_oi;
U1_oi=I1-R_oi; U2_oi=I2-R_oi;
S_oi=I12-R_oi-U1_oi-U2_oi;


out.i1=i1;
out.i2=i2;
out.i12=i12;
out.o123=o123;
out.icap_min=icap_min;
out.icap_oir=icap_oi;
out.I12=I12;
out.PRIDmin.R=R_min;
out.PRIDmin.S=S_min;
out.PRIDmin.U1=U1_min;
out.PRIDmin.U2=U2_min;
out.PRIDoir.R=R_oi;
out.PRIDoir.S=S_oi;
out.PRIDoir.U1=U1_oi;
out.PRIDoir.U2=U2_oi;
out.freq=outtmp.freq;


end