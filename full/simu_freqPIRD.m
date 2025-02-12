%% PARTIAL INFORMATION RATE DECOMPOSITION - coarse graining method k=1, M=3 sources
close all; clear; clc;

%% parameters 
M=4; % number of processes in the system
nsource=M-1;
Y=1; % target
X=setdiff(1:M,Y); % sources

coup=1; % coupling parameter
Fs=1; % sampling frequency
nfft=1000;
q=20; % number of lags for computing correlation functions

marksize=10; Dimfont=16;

range1=[0.04 0.15]; range2=[0.15 0.4];
nrange1=round((nfft*2/Fs)*range1);
if range2(2) < Fs/2
    nrange2=round((nfft*2/Fs)*range2)+[1 0];
else
    nrange2(1)=round((nfft*2/Fs)*range2(1))+1;
    nrange2(2)=round((nfft*2/Fs)*Fs/2);
end
band_all=[nrange1; nrange2];
nrange=size(band_all,1)+1;
band_name={'Time','B1','B2'};

%%% init
Red=zeros(1,nrange); Syn=Red;
Red_function_bottom_atom=Red;
Delta=Red;
Un=cell(1,nsource); I=Un;

delta='Δ';

%% theoretical simulation
par.Su=ones(1,M); % variance of innovation processes
par.poles{X(1)}=([0.8 0.3]); % source 1
par.poles{X(2)}=([0.8 0.3]); % source 2
par.poles{X(3)}=([0.9 0.1]); % source 3
par.coup=[X(1) Y 1 coup; ...
          X(1) X(2) 1 coup; ...
          X(3) Y 1 coup];

[Am,Su]=theoreticalVAR(M,par); % VAR parameters
Am=Am';

%% PIRD computation            
[Red_tmp,Syn_tmp,Delta_tmp,Un_tmp,MIR_tmp,...
Red_function_bottom_atom_tmp,i,r,f] = lrp_pird(Am,Su,q,Y,X,nfft,Fs,band_all);

for irange=1:nrange
    Red(1,irange)=Red_tmp{irange}; % redundancy
    Syn(1,irange)=Syn_tmp{irange}; % synergy
    Delta(1,irange)=Delta_tmp{irange}; % delta R-S
    Red_function_bottom_atom(1,irange)=Red_function_bottom_atom_tmp{irange}; % redundancy
    for im=1:nsource
        % PIRD
        Un{im}(1,irange) = Un_tmp{irange}(im); % unique X(im) --> Y
        I{im}(1,irange) = MIR_tmp{irange}(im); % MIR X(im) --> Y
    end

    % exemplary disp for fixed values of the parameters
    disp(['MMI redundancy; band ' band_name{irange}])
    disp(['Redundancy: ', num2str(Red(1,irange))])
    disp(['Synergy: ', num2str(Syn(1,irange))])
    disp([delta, ': ' num2str(Red(1,irange) ...
        - Syn(1,irange))])

    for im=1:nsource
        disp(['Unique X' num2str(im) '→ Y: ', num2str(Un{im}(1,irange))])
        disp(['MIR X' num2str(im) '→ Y: ', num2str(I{im}(1,irange))])
    end
    
end

%% plot 
% a figure is generated where the spectral profiles of MIR
% and the time-B1-B2 values of MIR, UNIQUE, RED, SYN are depicted as barplots

% figure('Color','w','WindowState','maximized');
figure('Color','w'); 
col_MIR{1}=[0 0 0]; col_MIR{2}=[0.5 0.5 0.5]; col_MIR{3}=[239/255 162/255 203/255];
col_red=[1 0 0]; col_syn=[0 0 1];

%%% MIR freq
subplot(2,2,1) % spectral profiles of the pairwise MIR
plot(f,r(8,:),'LineWidth',2,'Color',col_MIR{1}); hold on;
plot(f,r(9,:),'LineWidth',2,'Color',col_MIR{2}); 
plot(f,r(10,:),'LineWidth',2,'Color',col_MIR{3}); 
legend('i_{Y;X_1}','i_{Y;X_2}','i_{Y;X_3}'); legend box off
xlim([0 Fs/2]); 
xlabel('f [Hz]');
set(gca,'FontName','Times','FontSize',Dimfont);

%%% MIR time
subplot(2,2,2) % barplot: time-B1-B2 values of the pairwise MIR
n_MIR=3;
x_time=1:5:n_MIR*4+n_MIR;
x_freq{1}=2:5:n_MIR*4+n_MIR;
x_freq{2}=3:5:n_MIR*4+n_MIR;
k=1;
for im=1:nsource
    % time domain
    bar(x_time(k),I{im}(1,1),'EdgeColor',col_MIR{im},...
            'FaceColor',col_MIR{im}); hold on;
    % freq domain
    for irange=2:nrange
        bar(x_freq{irange-1}(k),I{im}(1,irange),'EdgeColor',col_MIR{im},...
        'FaceColor',col_MIR{im}); hold on;
    end
    k=k+1;
end
xticks([2 7 12]); xticklabels({'I_{Y;X_1}','I_{Y;X_2}','I_{Y;X_3}'});
xlim([0 x_freq{2}(end)+1])
set(gca,'FontName','Times','FontSize',Dimfont);

%%% MMI PIRD
subplot(2,2,[3 4]) % barplot: time-B1-B2 values of UN, RED & SYN
x_time=1:5:n_MIR*8+n_MIR;
x_freq{1}=2:5:n_MIR*8+n_MIR;
x_freq{2}=3:5:n_MIR*8+n_MIR;
k=1;
for im=1:nsource
    % time domain
    bar(x_time(k),Un{im}(1,1),'EdgeColor',col_MIR{im},...
            'FaceColor',col_MIR{im}); hold on;
    % freq domain
    for irange=2:nrange
        bar(x_freq{irange-1}(k),Un{im}(1,irange),'EdgeColor',col_MIR{im},...
        'FaceColor',col_MIR{im}); hold on;
    end
    k=k+1;
end
% redundancy
bar(x_time(k),Red(1,1),'EdgeColor',col_red,...
        'FaceColor',col_red); hold on;
for irange=2:nrange
    bar(x_freq{irange-1}(k),Red(1,irange),'EdgeColor',col_red,...
            'FaceColor',col_red); hold on;
end
% synergy
bar(x_time(k+1),Syn(1,1),'EdgeColor',col_syn,...
        'FaceColor',col_syn); hold on;
for irange=2:nrange
    bar(x_freq{irange-1}(k+1),Syn(1,irange),'EdgeColor',col_syn,...
            'FaceColor',col_syn); hold on;
end
xticks([2 7 12 17 22]); xticklabels({'U_{Y;X_1}','U_{Y;X_2}','U_{Y;X_3}','Red','Syn'});
xlim([0 x_freq{2}(end)+1])
set(gca,'FontName','Times','FontSize',Dimfont);



