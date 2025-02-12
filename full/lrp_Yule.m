%% Solution of Yule-Walker Equations for a VAR process (using discrete time Lyapunov eq)

function R = lrp_Yule(Am,Su,q)

%   Computes Information Dynamics analytically for a stationary mvar(p) process:
%   X_n=A_1*X_{n-1}+A_2*X_{t-n}+...+A_p*X_{n-p}+E_n
%
%   INPUT: 
%   Am  -   generalized connectivity matrix A=(A_1 A_2 ... A_p)
%   Su  -   covariance matrix for E_n
%   q   -   number of lags for which to compute correlations

%   OUTPUT:
%   R, QxQx(q+1) matrix of process covariances for lags k=0,1,...,q


Q = size(Am,1); %number of elements in the system
p=floor(size(Am,2)/Q); %number of lags in MVAR model

R=NaN*ones(Q,Q,q+1); % prepare covariance matrices, (:,:,1) is lag 0, (:,:,q+1) is lag q 

% Obtain F and Delta
Im=eye(Q*p);
Psi=[Am;Im(1:end-Q,:)];
Theta=zeros(p*Q,p*Q);
Theta(1:Q,1:Q)=Su(:,:);

% Obtain BigSigma solving the Lyapunov equation: BigSigma = A * BigSigma * A^T + Theta
BigSigma=dlyap(Psi,Theta);

% extract R(0),...,R(p-1)
for i=1:p
    R(:,:,i)=BigSigma(1:Q,Q*(i-1)+1:Q*i);
end

% Yule-Walker solution  for lags >= p
for k=p+1:q+1
    Rk=R(:,:,k-1:-1:k-p);
    Rm=[];
    for ki=1:p
        Rm=[Rm; Rk(:,:,ki)];
    end
    R(:,:,k)=Am*Rm;
end



