%% Computation of dynamic MIR and CMIR (linear Gaussian) for all pairs and blocks of series in a multivariate process
% --- given VAR parameters, computes MIR between each pair of processes, and conditional MIR given all other processes
% --- MIR and CMIR are computed also in the freq domain

% q: number of lags used to represent the past states of the processes
% Am, Su: VAR model parameters (theoretical or estimated with lrp_idVAR)
% nfft, Fs: number of points of freq axis and sampling freq

function out = lrp_MIRf_CMIRf(Am,Su,q,nfft,Fs)

Q=size(Su,1); % n. of processes

% initialize MIR and cMIR btw x and y
Ixy=zeros(Q,Q); Ixy_z=Ixy; 
Fxy=Ixy; Fxy_z=Ixy; 
% initialize spectral MIR and cMIR btw x and y
Ixyf=zeros(Q,Q,nfft); Ixy_zf=Ixyf; 
% initialize (spectral) MIR btw x and z
Ixz=zeros(Q,Q); Ixzf=zeros(Q,Q,nfft); 
% initialize (spectral) MIR btw blocks of processes
Ixyz=zeros(1,Q); Ixyzf=zeros(Q,nfft); 

for ix=1:Q
    for iy=ix+1:Q

        % compute MIR and MIRf
        iz=setdiff(1:Q,[ix iy]);
        retx_y=lrp_MIRf(Am,Su,q,ix,iy,nfft,Fs);
        Ixy(ix,iy)=retx_y.Ixy; % MIR btw processes ix and iy
        Ixyf(ix,iy,:)=retx_y.fxy; % MIRf btw processes ix and iy
        Fxy(ix,iy,:)=retx_y.Fxy; % integral of fxy
        f=retx_y.f; % freq axis

        % compute cMIR and cMIRf - btw processes ix and iy given processes iz 
        retx_z=lrp_MIRf(Am,Su,q,ix,iz,nfft,Fs); % MIR btw processes ix and iz
        retx_yz=lrp_MIRf(Am,Su,q,ix,[iy iz],nfft,Fs);  % MIR btw processes ix and [iy,iz]
        Ixy_z(ix,iy)=retx_yz.Ixy-retx_z.Ixy; % cMIR   
        Ixy_zf(ix,iy,:)=retx_yz.fxy-retx_z.fxy; % cMIRf  
        Fxy_z(ix,iy,:)=retx_yz.Fxy-retx_z.Fxy; % cMIR from the difference btw integrated values of cMIRf

    end

    % VAR model: Y; X1,X2,X3
    ioth=setdiff(1:Q,ix); % what's left (a block of processes)
    retxyz=lrp_MIRf(Am,Su,q,ix,ioth,nfft,Fs);  
    Ixyz(ix)=retxyz.Ixy; % MIR btw processes ix and ioth
    Ixyzf(ix,:)=retxyz.fxy; % MIRf btw processes ix and ioth

    % VAR model: Y; X(ix),X(iz) for all (ix,iz) in (1,2,...,M)
    if Q > 3 % this is run for more than 2 sources
        comb=nchoosek(1:length(ioth),2);
        for iz=1:length(ioth)
            isel=ioth(comb(iz,:));
            retxz=lrp_MIRf(Am,Su,q,ix,isel,nfft,Fs);  
            index=setdiff(1:Q,[ix isel]);
            Ixz(ix,index)=retxz.Ixy; % MIR
            Ixzf(ix,index,:)=retxz.fxy; % MIRf
        end
    end

end

%% OUTPUT
out.Ixy=Ixy; out.Ixyf=Ixyf; % MIR, MIRf (x,y)
out.Ixy_z=Ixy_z; out.Ixy_zf=Ixy_zf; % CMIR, CMIRf (x,y)

out.Fxy=Fxy; % MIR from integration
out.Fxy_z=Fxy_z; % CMIR from integration

if Q > 3
    out.Ixz=Ixz; out.Ixzf=Ixzf; % MIR, MIRf (x,z)
end
out.Ixyz=Ixyz; out.Ixyzf=Ixyzf; % MIR, MIRf (x,y,z)

out.f=f;

