%% MIR from linear regression
% performs linear regression of the present state of N target processes from the past states of M driver processes
% i: indexes of the first process, X
% j: indexes of the second process, Y
% q: number of lags used to represent the past states of the processes
% Am, Su: VAR model parameters (theoretical or estimated with lrp_idVAR)
% nfft, Fs: number of points of freq axis and sampling freq

function ret = lrp_MIRf(Am,Su,q,i,j,nfft,Fs)

M=length(i)+length(j);

%% TIME DOMAIN COMPUTATION
%%%% CE decomposition of MIR
% AR restricted models
[SigmaY,SigmaY_Y] = lrp_LinReg(Am,Su,q,j,j); % restricted model, Ypresent given Ypast
[SigmaX,SigmaX_X] = lrp_LinReg(Am,Su,q,i,i); % restricted model, Xpresent given Xpast
% full model
[SigmaXY,SigmaXY_XY,~,~,~,BXY] = lrp_LinReg(Am,Su,q,[i j],[i j]); % restricted model, XYpresent given XYpast

Hx=0.5*log(((2*pi*exp(1))^size(SigmaX,1))*det(SigmaX_X)); %entropy rate of X
Hy=0.5*log(((2*pi*exp(1))^size(SigmaY,1))*det(SigmaY_Y)); %entropy rate of Y
Hxy=0.5*log(((2*pi*exp(1))^size(SigmaXY,1))*det(SigmaXY_XY)); %entropy rate of [X Y]
Ixy=Hx+Hy-Hxy; % MIR - dec1

%%%% TE decomposition of MIR
% ARX restricted models
[~,SigmaY_XY] = lrp_LinReg(Am,Su,q,j,[i j]); % restricted model, Ypresent given Ypast and Xpast
[~,SigmaX_XY] = lrp_LinReg(Am,Su,q,i,[i j]); % restricted model, Xpresent given Ypast and Xpast

Txy=0.5*log(det(SigmaY_Y)/det(SigmaY_XY)); % TE from X to Y
Tyx=0.5*log(det(SigmaX_X)/det(SigmaX_XY)); % TE from Y to X
Ixoy=0.5*log(det(SigmaX_XY)*det(SigmaY_XY)/det(SigmaXY_XY)); % instantaneous TE
Ixy2=Txy+Tyx+Ixoy; % MIR - dec2

%% FREQUENCY DOMAIN COMPUTATION
% --- stability check ---
p = size(BXY,2)./M;
E=eye(M*p); AA=[BXY;E(1:end-M,:)];
lambda=eig(AA);
lambdamax=max(abs(lambda));
if lambdamax>=1 % The 2AR process is not stable 
    warning('The process is not stable');
end

[P,H,f] = lrp_VARspectra(BXY,SigmaXY_XY,nfft,Fs); %%% full model
P=abs(P); H=abs(H);
Px=P(1,1,:); % "i" is in position 1
Py=P(2:M,2:M,:); % "j" is in positions 2:M
Hx=H(1,1,:); Hy=H(2:M,2:M,:); 
fxy=zeros(1,nfft); %%% spectral MIR
fx_y=zeros(1,nfft); fy_x=zeros(1,nfft); %%% spectral TEs
for n=1:nfft
    dPx=abs(det(Px(:,:,n)));
    dPy=abs(det(Py(:,:,n)));
    dP=abs(det(P(:,:,n)));
    %%% spectral MIR
    fxy(n) = 0.5*log( (dPx .* dPy) ./ dP ); 
    %%% spectral TEs
    fx_y(n) = 0.5*log( dPy / abs(det(Hy(:,:,n)*SigmaXY_XY(2:M,2:M)*Hy(:,:,n)'))); % from x to y
    fy_x(n) = 0.5*log( dPx / abs(det(Hx(:,:,n)*SigmaXY_XY(1,1)*Hx(:,:,n)'))); % from y to x
end

%% VERIFY MIR
Fxy=sum(fxy)/nfft; % it should coincide with Ixy
Fx_y=sum(fx_y)/nfft; % it should coincide with Txy
Fy_x=sum(fy_x)/nfft; % it should coincide with Tyx

%% OUTPUT
ret.Ixy=Ixy; % MIR
ret.fxy=fxy; % MIRf
ret.Fxy=Fxy; % integral of fxy

ret.Txy=Txy; % TE x-->y
ret.fx_y=fx_y; % freq. domain TE
ret.Fx_y=Fx_y; % integral of fx_y

ret.Tyx=Tyx; % TE y-->x
ret.fy_x=fy_x; % freq. domain TE
ret.Fy_x=Fy_x; % integral of fy_x

ret.f=f; % freq
