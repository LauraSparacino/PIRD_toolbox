%% PIRD and PID frameworks applied to climate data
clear; close all; clc;

% addpath([pwd '\functions\']);
% addpath([pwd '\data\']);

%% parameters
nfft=1000; %number of points on frequency axis 
fs=1; % Sampling frequency (1 month)
pfilter=0.94; % filter parameter for detrending
q=20; % for computing correlation functions
pmax=14; % maximum scanned model order (Akaike Information Criterion - AIC)
alpha=0.05; % threshold for surrogate data analysis
Ns=100; % number of surrogates
M=3; % number of network units
nsource=M-1; % number of sources

% plot parameters:
diff=0.3;
msize=4;
msize_ast=8;
pos{1}=[1 5]; pos{2}=[2 6]; pos{3}=[3 7]; 

%% load data
load("climate_data.mat");
[Nsamples,indend]=size(clima);

%% sources are scalar 
isource=nan(1,nsource);

itarget=7; % SOI - target series
% select the sources among the 13 series: NINO34, TSA, PDO, NTA (11-8-10-12)
isource_tmp=[11 8;...
            11 10;...
            11 12;...
            8 10;...
            8 12;...
            10 12];
ncomb=size(isource_tmp,1);

%figure('WindowState','maximized','Color','w'); 
figure('Color','w'); 

for ic=1:ncomb
    % for each combination (for each triplet), build Ns surrogates (random
    % shuffling) to destroy temporal correlations but maintain zero-lag
    % interactions (i.e., PIRD varies but PID remains the same)
   
    disp(['Combination: ', num2str(ic)])
    
    isource=isource_tmp(ic,:);
    idata=[itarget isource];
    data_o=clima(:,idata); % original data

    names=indici_clima(idata);
    name_target=names{1};
    name_s1=names{2}; name_s2=names{3};

    % pre-processing: filtering and normalization
    data_f=data_o; data=data_o;
    for ii=1:M
        data_f(:,ii)=AR_filter(data_o,ii,pfilter); % AR highpass filtered series
        data(:,ii)=zscore(data_f(:,ii)); % zero-mean series with unitary variance
    end

    % cycle on surrogates (the first cycle is on original data)
    for ns=1:Ns+1
        
        if ns==1
            data_surr=data;
        else
            rp=randperm(Nsamples);
            data_surr=zeros(Nsamples,M);
            for k = 1:Nsamples
                data_surr(k,:)=data(rp(k),:);
            end
        end

        %% data analysis --- PIRD 
        % model order selection - model identification
        p = mos_idMVAR(data_surr',pmax,0);
        out=lrp_idVAR(data_surr,1:M,1:M,1:p);
        Am=out.eA; Am=Am';
        Su=out.es2u;
        Y=1; X=2:M;
        
        [R_pird_tmp,S_pird_tmp,Delta_pird_tmp,U_pird_tmp,I_pird_tmp,...
            ~,i_pird,r_pird,f] = lrp_pird(Am,Su,q,Y,X,nfft,fs);
        R_pird(ns)=R_pird_tmp{1}; % time domain values
        S_pird(ns)=S_pird_tmp{1}; 
        U1_pird(ns)=U_pird_tmp{1}(1);
        U2_pird(ns)=U_pird_tmp{1}(2); 
        I1_pird(ns)=I_pird_tmp{1}(1);
        I2_pird(ns)=I_pird_tmp{1}(2);
        I12_pird(ns)=I_pird_tmp{1}(3);

        %% data analysis --- PID 
        %%% from time series: static model
        out=lrp_LinRegStatic(data_surr,Y,X); % full model
        I12_pid(ns)=-0.5*log(1-out.erho2);
        out=lrp_LinRegStatic(data_surr,Y,X(1)); % restricted model: X1
        I1_pid(ns)=-0.5*log(1-out.erho2);
        out=lrp_LinRegStatic(data_surr,Y,X(2)); % restricted model: X2
        I2_pid(ns)=-0.5*log(1-out.erho2);
        R_pid(ns)=min([I1_pid(ns) I2_pid(ns)]); % zero-lag MMI redundancy
        U1_pid(ns)=I1_pid(ns)-R_pid(ns);
        U2_pid(ns)=I2_pid(ns)-R_pid(ns);
        S_pid(ns)=I12_pid(ns)-U1_pid(ns)-U2_pid(ns)-R_pid(ns);
        
    end
       
    subplot(1,5,1)
    plot(ic+diff,I12_pid(1),'*','Markersize',msize_ast,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on;
    plot(ic+diff,I12_pird(1),'*','MarkerFaceColor','k','MarkerEdgeColor','k',...
        'Markersize',msize_ast); hold on;
    plot(ic,I12_pird(2:end),'o','MarkerFaceColor','k','MarkerEdgeColor','k',...
        'Markersize',msize); hold on;
    xlim([0 7])
    xticks(1:ncomb); xticklabels({'1)','2)','3)','4)','5)','6)'});
    title('joint MIR');
    pbaspect([1 1 1])
    subplot(1,5,2)
    plot(ic+diff,U1_pid(1),'*','Markersize',msize_ast,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on;
    plot(ic+diff,U1_pird(1),'*','Markersize',msize_ast,'MarkerFaceColor',[0.4660 0.6740 0.1880],...
        'MarkerEdgeColor',[0.4660 0.6740 0.1880]); hold on;
    plot(ic,U1_pird(2:end),'o','Markersize',msize,'MarkerFaceColor',[0.4660 0.6740 0.1880],...
        'MarkerEdgeColor',[0.4660 0.6740 0.1880]); hold on;
    title('Unique 1');
    xlim([0 7])
    xticks(1:ncomb); xticklabels({'1)','2)','3)','4)','5)','6)'});
    pbaspect([1 1 1])
    subplot(1,5,3)
    plot(ic+diff,U2_pid(1),'*','Markersize',msize_ast,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on;
    plot(ic+diff,U2_pird(1),'*','MarkerFaceColor',[0.85 0.325 0.098],'MarkerEdgeColor',[0.85 0.325 0.098],...
        'Markersize',msize_ast); hold on;
    plot(ic,U2_pird(2:end),'o','MarkerFaceColor',[0.85 0.325 0.098],'MarkerEdgeColor',[0.85 0.325 0.098],...
        'Markersize',msize); hold on;
    title('Unique 2');
    xlim([0 7])
    xticks(1:ncomb); xticklabels({'1)','2)','3)','4)','5)','6)'});
    pbaspect([1 1 1])
    subplot(1,5,4)
    plot(ic+diff,S_pid(1),'*','Markersize',msize_ast,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on;
    plot(ic+diff,S_pird(1),'*','MarkerFaceColor','b','MarkerEdgeColor','b',...
        'Markersize',msize_ast); hold on;
    plot(ic,S_pird(2:end),'o','MarkerFaceColor','b','MarkerEdgeColor','b',...
        'Markersize',msize); hold on;
    title('Synergy');
    xlim([0 7])
    xticks(1:ncomb); xticklabels({'1)','2)','3)','4)','5)','6)'});
    pbaspect([1 1 1])
    subplot(1,5,5)
    plot(ic+diff,R_pid(1),'*','Markersize',msize_ast,'MarkerEdgeColor',[0.5 0.5 0.5]); hold on;
    plot(ic+diff,R_pird(1),'*','MarkerFaceColor','r','MarkerEdgeColor','r',...
        'Markersize',msize_ast); hold on;
    plot(ic,R_pird(2:end),'o','MarkerFaceColor','r','MarkerEdgeColor','r',...
        'Markersize',msize); hold on;
    title('Redundancy');
    xlim([0 7])
    xticks(1:ncomb); xticklabels({'1)','2)','3)','4)','5)','6)'});
    pbaspect([1 1 1])
end

% rmpath([pwd '\functions\']);
% rmpath([pwd '\data\']);

