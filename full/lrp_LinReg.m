%% Linear regression of random processes
% performs linear regression of the present state of N target processes from the past states of M driver processes
% j: indexes of target processes (processes to be predicted)
% i: indexes of driver processes (predictors)
% q: number of lags used to represent the past states of the processes
% Am, Su: VAR model parameters (theoretical or estimated with lrp_idVAR)

function [SigmaY,SigmaY_X,SigmaS,SigmaX,SigmaYX,B,R] = lrp_LinReg(Am,Su,q,j,i)

N=length(j);
M=length(i);

R = lrp_Yule(Am,Su,q); % Lambda(0),...,Lambda(q) - dim: M*M*(q+1)

SigmaS=R(:,:,1); % covariance of all processes
SigmaY=SigmaS(j,j); % covariance of target processes

SigmaX=NaN*ones(M*q,M*q); %covariance of the past of the driver processes
PR=reshape(R(i,i,1:q),[M,M*q]);
SigmaX(1:M,:)=PR;
SigmaX(:,1:M)=PR';
for k=1:q-1
    PR(:,end-M+1:end)=[];
    SigmaX(k*M+1:k*M+M,end-size(PR,2)+1:end)=PR;
    SigmaX(end-size(PR,2)+1:end,k*M+1:k*M+M)=PR';
end

SigmaYX=NaN*ones(N,M*q); % cross-covariance between the present of the target processes and the past of the driver processes
for k=1:q
    SigmaYX(:,M*(k-1)+1:M*(k-1)+M)=R(j,i,k+1);
end

B=SigmaYX/SigmaX; % B=SigmaYX*inv(SigmaX);

SigmaY_X=SigmaY-B*SigmaX*B'; % partial covariance of the target given the driver processes


