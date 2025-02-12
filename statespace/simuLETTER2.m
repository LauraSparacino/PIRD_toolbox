%% test PIRD of random processes vs PID of random variables
% PIRD/PID with one target process Y and two source processes X1,X2

clear; close all; clc;

%%% SIMULATION setting
M=3; p=1; Mv=[1 1 1];
d=0:0.02:1;
N=M-1;

Fs=1; nfft=512;
it=3; %index of target
is=setdiff(1:M,it); %index of sources

%% analysis cycle

for i=1:length(d)
    clc; disp(['simu ' int2str(i) ' of ' int2str(length(d))])
    
    c1=d(i);
    c2=d(i);
    a1=0.8*(1-d(i)/max(d)); a2=0.8*(1-d(i)/max(d));
    b1=0.1; b2=b1;
    Su=eye(M);
    
    Ak(:,:,1)=[a1 b1 0; b2 a2 0; c1 c2 0];
    Am=[]; %group coefficient blocks Ak in a single matrix Am
    for kk=1:p
        Am=[Am Ak(:,:,kk)];
    end
    % stability check
    E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
    if lambdamax>=1, error('The simulated VAR process is not stable'); end
    
    
    %%% analysis - PID
    R = lrp_Yule(Am,Su,p); % Lambda(0),...,Lambda(q) - dim: M*M*(q+1)
    SigmaS=R(:,:,1); % lsg-zero covariance of all processes
    SigmaY=SigmaS(it,it);
    SigmaX=SigmaS(is,is);
    SigmaX_Y=SigmaS(is,it);
    SigmaW = SigmaY - SigmaX_Y' * inv(SigmaX') * SigmaX_Y;
    %B=SigmaX_Y'/SigmaX; % B=SigmaYX*inv(SigmaX);
    %SigmaY_X=SigmaY-B*SigmaX*B'; % partial covariance of the target given the driver processes
    SigmaW1 = SigmaY - SigmaS(is(1),it)^2 / SigmaS(is(1),is(1));
    SigmaW2 = SigmaY - SigmaS(is(2),it)^2 / SigmaS(is(2),is(2));
    I_Y_X(i)=0.5*log(SigmaY/SigmaW);
    I_Y_X1(i)=0.5*log(SigmaY/SigmaW1);
    I_Y_X2(i)=0.5*log(SigmaY/SigmaW2);
    R_pid(i)=min([I_Y_X1(i) I_Y_X2(i)]); % zero-lag MMI redundancy
    U1_pid(i)=I_Y_X1(i)-R_pid(i);
    U2_pid(i)=I_Y_X2(i)-R_pid(i);
    S_pid(i)=I_Y_X(i)-U1_pid(i)-U2_pid(i)-R_pid(i);

    
    %%% analysis - PIRD
    [A,C,K,~] = oir_ar2iss(Am); %%% equivalent ISS model
    out=pird2(A,C,K,Su,Mv,it,Fs,nfft); % PIRD with two sources  
    I12(i)=out.I12;
    R_pird(i)=out.PRIDmin.R;
    U1_pird(i)=out.PRIDmin.U1;
    U2_pird(i)=out.PRIDmin.U2;
    S_pird(i)=out.PRIDmin.S;
end



%% plots
col1=[0 0 0]/255; % black
col1b=[96 96 96]/255; % dark grey
col1c=[192 192 192]/255; % light grey
col2=[192 64 64]/255; % dark red
col2b=[192 64 160]/255; % viola
col2c=[192 160 64]/255; % miele
col3=[64 128 192]/255; % dark azure
col3b=[64 64 192]/255; % blu
col3c=[64 192 192]/255; % celeste
col4=[64 160 64]/255; % dark green
col5=[224 224 32]/255; % yellow
lw=1.2;

figure(1)
subplot(1,2,1);
plot(d,I_Y_X,'color',col1b,'linewidth',lw); hold on;
plot(d,R_pid,'color',col2,'linewidth',lw);
plot(d,S_pid,'color',col4,'linewidth',lw,'linestyle','--');
plot(d,U1_pid,'color',col3,'linewidth',lw);
plot(d,U2_pid,'color',col3b,'linewidth',lw,'linestyle','--');
xlim([min(d) max(d)]);
legend('I(Y_n;X_{n})','R(Y_n;X_{n})','S(Y_n;X_{n})','U(Y_n;X_{1,n})','U(Y_n;X_{2,n})')
title('Zero-lag Partial Information Decomposition')
xlabel('d')

subplot(1,2,2);
plot(d,I12,'color',col1b,'linewidth',lw); hold on;
plot(d,R_pird,'color',col2,'linewidth',lw);
plot(d,S_pird,'color',col4,'linewidth',lw,'linestyle','--');
plot(d,U1_pird,'color',col3,'linewidth',lw);
plot(d,U2_pird,'color',col3b,'linewidth',lw,'linestyle','--');
xlim([min(d) max(d)]);
legend('I_{Y;X}','R_{Y;X}','S_{Y;X}','U_{Y;X_1}','U_{Y;X_2}')
title('Partial Information Rate Decomposition')
xlabel('d')

