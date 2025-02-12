%% PIRD decomposition using coarse graining approach (k=1) in case of N=2 and N=3 source processes

% Y: target
% X: source vector process
% band_all: (b x 2) matrix with b=number of freq. ranges
% natoms: number of atoms in the lattice structure
% q: number of lags used to represent the past states of the processes
% Am, Su: VAR model parameters (theoretical or estimated with lrp_idVAR)
% nfft, Fs: number of points of freq axis and sampling freq
% do_freq: 'y' or 'n' for spectral analysis in given frequency bands or not

function [Red,Syn,Delta,Un,MIR,...
    Red_function_bottom_atom,i,r,f] = lrp_pird(Am,Su,q,Y,X,nfft,Fs,band_all)

% check input: number of atoms of the lattice
if length(X)==2; natoms=4;
elseif length(X)==3; natoms=18;
else; error('Sorry, no lattice for more than 3 sources...')
end
% check input: time/frequency domain analysis
if nargin<8; band_all=nan; % if not specified, band_all is nan and only time domain analysis is performed
elseif nargin<7; Fs=1; 
elseif nargin<6; nfft=1000;
end

if ~isnan(band_all)
    b=size(band_all,1)+1; % number of freq bands 
    do_freq='y';
else
    b=1; % compute only time domain PIRD measures
    do_freq='n';
end

% init
r=nan*ones(natoms,nfft); % spectral atom redundancy 
i=nan*ones(natoms,nfft); % spectral atom information 
Red=cell(1,b); Syn=Red; Un=Red; Delta=Red; MIR=Red; % PIRD 
% redundancy functions associated to the bottom atom of the lattice:
Red_function_bottom_atom=cell(1,b);

%% computation --- MIR between blocks of processes
% this computes dynamic MIR and CMIR (linear Gaussian) for all pairs and blocks of series 
% computation is performed in the time and spectral domains
out=lrp_MIRf_CMIRf(Am,Su,q,nfft,Fs); 

% MIR between pairs of processes (Y;X(i))
Ixyf=out.Ixyf; % spectral profile
% MIR between Y and the block of all the other sources
Ixyzf=out.Ixyzf; % spectral profiles
f=out.f; % frequency axis

%% lattice with 4 atoms, 2 sources
if natoms==4 % N=2 sources
    
    for nf=1:nfft
        
        %%% SPECTRAL ATOM REDUNDANCY
        r(1,nf)=min(Ixyf(Y,X(1),nf),Ixyf(Y,X(2),nf)); % r{1}{2}=min(i(Y;X1),i(Y;X2))=icap{1}{2}
        r(2,nf)=Ixyf(Y,X(1),nf); % r{1}=i(Y;X1)=icap{1}
        r(3,nf)=Ixyf(Y,X(2),nf); % r{2}=i(Y;X2)=icap{2}
        r(4,nf)=Ixyzf(Y,nf); % r{12}=i(Y;X1,X2)=icap{12}
        %%% SPECTRAL ATOM INFORMATION
        i(1,nf)=r(1,nf); % redundancy
        i(2,nf)=r(2,nf)-r(1,nf); % U1
        i(3,nf)=r(3,nf)-r(1,nf); % U2
        i(4,nf)=r(4,nf)-(i(1,nf)+i(2,nf)+i(3,nf)); % synergy
        
    end
    
    % time domain 
    Red{1}=sum(i(1,:))/nfft;
    Un{1}(1)=sum(i(2,:))/nfft;
    Un{1}(2)=sum(i(3,:))/nfft;
    Syn{1}=sum(i(4,:))/nfft;
    Delta{1}=Red{1}-Syn{1};
    MIR{1}(1)=sum(r(2,:))/nfft; % Y, X1
    MIR{1}(2)=sum(r(3,:))/nfft; % Y, X2
    MIR{1}(3)=sum(r(4,:))/nfft; % Y, [X1,X2]
    Red_function_bottom_atom{1}=Red{1};
    % frequency domain
    if do_freq~='n'
        for irange=2:b
            %%% COARSE-GRAINED PIRD 
            Red{irange}=sum(i(1,band_all(irange-1,1):band_all(irange-1,2)))/nfft; 
            Un{irange}(1)=sum(i(2,band_all(irange-1,1):band_all(irange-1,2)))/nfft; 
            Un{irange}(2)=sum(i(3,band_all(irange-1,1):band_all(irange-1,2)))/nfft; 
            Syn{irange}=sum(i(4,band_all(irange-1,1):band_all(irange-1,2)))/nfft; 
            Delta{irange}=Red{irange}-Syn{irange};
            %%% MUTUAL INFORMATION 
            MIR{irange}(1)=sum(r(2,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % Y, X1
            MIR{irange}(2)=sum(r(3,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % Y, X2
            MIR{irange}(3)=sum(r(4,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % Y, [X1,X2]
            %MIR_test(1)=Un(1)+Red; % Y, X1
            %MIR_test(2)=Un(2)+Red; % Y, X2
            Red_function_bottom_atom{irange}=Red{irange};
        end
    end
    
%% lattice with 18 atoms, 3 sources
elseif natoms==18 % N=3 sources
    
    % MIR between Y and blocks of sources; e.g. I(Y;X1,X2) is in position Y4 
    Ixzf=out.Ixzf; % spectral profile
    
    for nf=1:nfft

        %%% --- SPECTRAL ATOM REDUNDANCY
        % 1st row
        r(1,nf)=min([Ixyf(Y,X(1),nf) Ixyf(Y,X(2),nf) Ixyf(Y,X(3),nf)]); lattice{1}='{1}{2}{3}'; %icap{1}{2}{3}
        % 2nd row
        r(2,nf)=min([Ixyf(Y,X(1),nf) Ixyf(Y,X(2),nf)]); lattice{2}='{1}{2}'; %icap{1}{2}
        r(3,nf)=min([Ixyf(Y,X(1),nf) Ixyf(Y,X(3),nf)]); lattice{3}='{1}{3}'; %icap{1}{3}
        r(4,nf)=min([Ixyf(Y,X(2),nf) Ixyf(Y,X(3),nf)]); lattice{4}='{2}{3}'; %icap{2}{3}
        % 3rd row
        r(5,nf)=min([Ixyf(Y,X(1),nf) Ixzf(Y,X(1),nf)]); lattice{5}='{1}{23}'; %icap{1}{23}
        r(6,nf)=min([Ixyf(Y,X(2),nf) Ixzf(Y,X(2),nf)]); lattice{6}='{2}{13}'; %icap{2}{13}
        r(7,nf)=min([Ixyf(Y,X(3),nf) Ixzf(Y,X(3),nf)]); lattice{7}='{3}{12}'; %icap{3}{12}
        % 4th row
        r(8,nf)=Ixyf(Y,X(1),nf); lattice{8}='{1}'; %icap{1}
        r(9,nf)=Ixyf(Y,X(2),nf); lattice{9}='{2}'; %icap{2}
        r(10,nf)=Ixyf(Y,X(3),nf); lattice{10}='{3}'; %icap{3}
        % 5th row
        r(11,nf)=min([Ixzf(Y,X(1),nf) Ixzf(Y,X(2),nf) Ixzf(Y,X(3),nf)]); lattice{11}='{12}{13}{23}'; %icap{12}{13}{23}
        % 6th row
        r(12,nf)=min([Ixzf(Y,X(3),nf) Ixzf(Y,X(2),nf)]); lattice{12}='{12}{13}'; %icap{12}{13}
        r(13,nf)=min([Ixzf(Y,X(3),nf) Ixzf(Y,X(1),nf)]); lattice{13}='{12}{23}'; %icap{12}{23}
        r(14,nf)=min([Ixzf(Y,X(2),nf) Ixzf(Y,X(1),nf)]); lattice{14}='{13}{23}'; %icap{13}{23}
        % 7th row
        r(15,nf)=Ixzf(Y,X(3),nf); lattice{15}='{12}'; %icap{12}
        r(16,nf)=Ixzf(Y,X(2),nf); lattice{16}='{13}'; %icap{13}
        r(17,nf)=Ixzf(Y,X(2),nf); lattice{17}='{23}'; %icap{23}
        % 8th row
        r(18,nf)=Ixyzf(Y,nf); lattice{18}='{123}'; %icap{123}
        
        %%% --- SPECTRAL ATOM INFORMATION    
        % 1st row
        i(1,nf)=r(1,nf);
        % 2nd row
        i(2,nf)=r(2,nf)-i(1,nf);
        i(3,nf)=r(3,nf)-i(1,nf);
        i(4,nf)=r(4,nf)-i(1,nf);
        % 3rd row
        i(5,nf)=r(5,nf)-sum(i([1 2 3],nf));
        i(6,nf)=r(6,nf)-sum(i([1 2 4],nf));
        i(7,nf)=r(7,nf)-sum(i([1 3 4],nf));
        % 4th row
        i(8,nf)=r(8,nf)-sum(i([1 2 3 5],nf));
        i(9,nf)=r(9,nf)-sum(i([1 2 4 6],nf));
        i(10,nf)=r(10,nf)-sum(i([1 3 4 7],nf));
        % 5th row
        i(11,nf)=r(11,nf)-sum(i([1 2 3 4 5 6 7],nf));
        % 6th row
        i(12,nf)=r(12,nf)-sum(i([1 2 3 4 5 6 7 8 11],nf));
        i(13,nf)=r(13,nf)-sum(i([1 2 3 4 5 6 7 9 11],nf));
        i(14,nf)=r(14,nf)-sum(i([1 2 3 4 5 6 7 10 11],nf));
        % 7th row
        i(15,nf)=r(15,nf)-sum(i([1 2 3 4 5 6 7 8 9 11 12 13],nf));
        i(16,nf)=r(16,nf)-sum(i([1 2 3 4 5 6 7 8 10 11 12 14],nf));
        i(17,nf)=r(17,nf)-sum(i([1 2 3 4 5 6 7 9 10 11 13 14],nf));
        % 8th row
        i(18,nf)=r(18,nf)-sum(i(1:17,nf));
        
    end
            
    % time domain     
    Red_function_bottom_atom{1}=sum(r(1,:))/nfft;
    %%% COARSE-GRAINED PIRD
    Red{1}=sum(i(1,:))/nfft + sum(i(2,:))/nfft + sum(i(3,:))/nfft + sum(i(4,:))/nfft; % 1st and 2nd rows
    Un{1}(1)=sum(i(5,:))/nfft + sum(i(8,:))/nfft; % left 2 atoms of the 1st column (3rd and 4th rows)
    Un{1}(2)=sum(i(6,:))/nfft + sum(i(9,:))/nfft; % middle 2 atoms of the 2nd column (3rd and 4th rows)
    Un{1}(3)=sum(i(7,:))/nfft + sum(i(10,:))/nfft; % right 2 atoms of the 3rd column (3rd and 4th rows)
    Syn{1}=sum(i(11,:))/nfft + sum(i(12,:))/nfft + sum(i(13,:))/nfft + sum(i(14,:))/nfft + ...
        sum(i(15,:))/nfft + sum(i(16,:))/nfft + sum(i(17,:))/nfft + sum(i(18,:))/nfft; % all the remaining atoms
    Delta{1}=Red{1}-Syn{1};
    %%% MUTUAL INFORMATION 
    MIR{1}(1)=sum(r(8,:))/nfft; % Y, X1
    MIR{1}(2)=sum(r(9,:))/nfft; % Y, X2
    MIR{1}(3)=sum(r(10,:))/nfft; % Y, X3

    % frequency domain
    if do_freq~='n'
        for irange=2:b
            Red_function_bottom_atom{irange}=sum(r(1,band_all(irange-1,1):band_all(irange-1,2)))/nfft;
            %%% COARSE-GRAINED PIRD 
            Red{irange}=sum(i(1,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(2,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(3,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(4,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % 1st and 2nd rows
            Un{irange}(1)=sum(i(5,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(8,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % left 2 atoms of the 1st column (3rd and 4th rows)
            Un{irange}(2)=sum(i(6,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(9,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % middle 2 atoms of the 2nd column (3rd and 4th rows)
            Un{irange}(3)=sum(i(7,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(10,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % right 2 atoms of the 3rd column (3rd and 4th rows)
            Syn{irange}=sum(i(11,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(12,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(13,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(14,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(15,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(16,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(17,band_all(irange-1,1):band_all(irange-1,2)))/nfft + ...
                sum(i(18,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % all the remaining atoms
            Delta{irange}=Red{irange}-Syn{irange};
            %%% MUTUAL INFORMATION 
            MIR{irange}(1)=sum(r(8,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % Y, X1
            MIR{irange}(2)=sum(r(9,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % Y, X2
            MIR{irange}(3)=sum(r(10,band_all(irange-1,1):band_all(irange-1,2)))/nfft; % Y, X3
        end    
    end
    
else
    error('Check the number of atoms of the lattice and the number of sources.')
end

end