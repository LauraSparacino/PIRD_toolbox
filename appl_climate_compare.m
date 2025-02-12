%% PIRD and PID frameworks applied to climate data
%% COMPARISON of full framework based on submodels with covariance matrix (Laura) and state-space framework implemented only for 2 sources (Luca)
clear; close all; clc;

addpath([pwd '\full\']);
addpath([pwd '\statespace\']);

% parameters
itarget=7; % target series 7 = SOI
isource=[11 12];% select the sources among the 13 series: NINO34, TSA, PDO, NTA (11-8-10-12)
idata=[itarget isource];
M=length(idata); % nr of processes
pmax=14; % maximum scanned model order (Akaike Information Criterion - AIC)
q=20; % for computing correlation functions
nfft=1000; %number of points on frequency axis 
fs=1; % Sampling frequency (1 month)


%% open data, select series, preprocess
load("climate_data.mat");
[Nsamples,indend]=size(clima);
data_o=clima(:,idata); % original data

pfilter=0.94; % filter parameter for detrending
data_f=data_o; data=data_o;
for ii=1:M
    data_f(:,ii)=AR_filter(data_o,ii,pfilter); % AR highpass filtered series
    data(:,ii)=zscore(data_f(:,ii)); % zero-mean series with unitary variance
end

names=indici_clima(idata);
name_target=names{1}; name_s1=names{2}; name_s2=names{3};


%% data analysis --- PIRD 
% model order selection - model identification
p = mos_idMVAR(data',pmax,0);
out=lrp_idVAR(data,1:M,1:M,1:p);
Am=out.eA';
Su=out.es2u;

% PIRD (Laura)
iY=1; iX=2:M;      
[Red_r,Syn_r,Delta_r,Un_r,I_r,i_r,r_r,f] = lrp_pird(Am,Su,q,iY,iX,nfft,fs);
% time domain values
R_pird=Red_r{1}; 
S_pird=Syn_r{1}; 
U1_pird=Un_r{1}(1);
U2_pird=Un_r{1}(2); 
I1_pird=I_r{1}(1);
I2_pird=I_r{1}(2);
I12_pird=I_r{1}(3);

%%% state space PIRD (Luca)
it=1; %index of target
is=setdiff(1:M,it); %index of sources
Mv=[1 1 1];
[A,C,K,~] = oir_ar2iss(Am); %%% equivalent ISS model
out=pird2(A,C,K,Su,Mv,it,fs,nfft); % PIRD with two sources  
I12ss=out.I12;
R_sspird=out.PRIDmin.R;
U1_sspird=out.PRIDmin.U1;
U2_sspird=out.PRIDmin.U2;
S_sspird=out.PRIDmin.S;

%% data analysis --- PID 
%%% from time series: static model (Laura)
out=lrp_LinRegStatic(data,iY,iX); % full model
I12_pid=-0.5*log(1-out.erho2);
out=lrp_LinRegStatic(data,iY,iX(1)); % restricted model: X1
I1_pid=-0.5*log(1-out.erho2);
out=lrp_LinRegStatic(data,iY,iX(2)); % restricted model: X2
I2_pid=-0.5*log(1-out.erho2);
R_pid=min([I1_pid I2_pid]); % zero-lag MMI redundancy
U1_pid=I1_pid-R_pid;
U2_pid=I2_pid-R_pid;
S_pid=I12_pid-U1_pid-U2_pid-R_pid;


%%% PID based on Yule-Walker (Luca)
R = lrp_Yule(Am,Su,p); % Lambda(0),...,Lambda(q) - dim: M*M*(q+1)
SigmaS=R(:,:,1); % lag-zero covariance of all processes
SigmaY=SigmaS(it,it);
SigmaX=SigmaS(is,is);
SigmaX_Y=SigmaS(is,it);
SigmaW = SigmaY - SigmaX_Y' * inv(SigmaX') * SigmaX_Y;
%B=SigmaX_Y'/SigmaX; % B=SigmaYX*inv(SigmaX);
%SigmaY_X=SigmaY-B*SigmaX*B'; % partial covariance of the target given the driver processes
SigmaW1 = SigmaY - SigmaS(is(1),it)^2 / SigmaS(is(1),is(1));
SigmaW2 = SigmaY - SigmaS(is(2),it)^2 / SigmaS(is(2),is(2));
I_Y_X=0.5*log(SigmaY/SigmaW);
I_Y_X1=0.5*log(SigmaY/SigmaW1);
I_Y_X2=0.5*log(SigmaY/SigmaW2);
R_ywpid=min([I_Y_X1 I_Y_X2]); % zero-lag MMI redundancy
U1_ywpid=I_Y_X1-R_ywpid;
U2_ywpid=I_Y_X2-R_ywpid;
S_ywpid=I_Y_X-U1_ywpid-U2_ywpid-R_ywpid;

%% display
disp('Comparison PIRD results')
disp('submodels from cov matrix at lag q (Laura):')
disp(['total MIR:' num2str(I12_pird)]);
disp(['U1=' num2str(U1_pird) ', U2=' num2str(U2_pird) ',Red=' num2str(R_pird) ',Syn=' num2str(S_pird)])
disp('submodels from cstate space method (Luca):')
disp(['total MIR:' num2str(I12ss)]);
disp(['U1=' num2str(U1_sspird) ', U2=' num2str(U2_sspird) ',Red=' num2str(R_sspird) ',Syn=' num2str(S_sspird)])

disp('Comparison zero-lag PID results')
disp('identifying several static models (Laura):')
disp(['total MI:' num2str(I12_pid)]);
disp(['U1=' num2str(U1_pid) ', U2=' num2str(U2_pid) ',Red=' num2str(R_pid) ',Syn=' num2str(S_pid)])
disp('single model, solving YW eqs (Luca):')
disp(['total MI:' num2str(I_Y_X)]);
disp(['U1=' num2str(U1_ywpid) ', U2=' num2str(U2_ywpid) ',Red=' num2str(R_ywpid) ',Syn=' num2str(S_ywpid)])


rmpath([pwd '\full\']);
rmpath([pwd '\statespace\']);

