%% VAR Spectral matrices (cross-PSD, TF)

%%% inputs: 
% Am=[A(1)...A(p)]: M*pM matrix of the MVAR model coefficients (strictly causal model)
% Su: M*M covariance matrix of the input noises
% N= number of points for calculation of the spectral functions (nfft)
% fs= sampling frequency

%%% outputs:
% H= Tranfer Function Matrix 
% S= Spectral Matrix 
% f= frequency vector

function [S,H,f] = lrp_VARspectra(Am,Su,N,fs)

M= size(Am,1); % Am has dim M*pM
p = size(Am,2)/M; % p is the order of the MVAR model

if nargin<2, Su = eye(M,M); end % if not specified, we assume uncorrelated noises with unit variance as inputs 
if nargin<3, N = 512; end
if nargin<4, fs= 1; end   
if all(size(N)==1)	 %if N is scalar
    f = (0:N-1)*(fs/(2*N)); % frequency axis
else            % if N is a vector, we assume that it is the vector of the frequencies
    f = N; N = length(N);
end

% s = exp(1i*2*pi*f/fs); % vector of complex exponentials
z = 1i*2*pi/fs;

% Initializations: spectral matrices have M rows, M columns and are calculated at each of the N frequencies
H=zeros(M,M,N); % Transfer Matrix
S=zeros(M,M,N); % Spectral Matrix

A = [eye(M) -Am]; % matrix from which M*M blocks are selected to calculate spectral functions
% Su=diag(diag(Su));

%% computation of spectral functions
for n=1:N % at each frequency
    
        %%% Coefficient matrix in the frequency domain
        As = zeros(M,M); % matrix As(z)=I-sum(A(k))
        for k = 1:p+1
            As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n)); 
            % indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix A (A(1) is in the second block, and so on)
        end
        
        %%% Transfer matrix 
        H(:,:,n)  = inv(As);
        
        %%% Spectral matrix 
        S(:,:,n)  = H(:,:,n)*Su*H(:,:,n)'; % ' stands for Hermitian transpose
       
end



