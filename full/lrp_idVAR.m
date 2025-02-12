%% GENERAL LINEAR REGRESSION MODEL
% given data matrix and indexes of predicted and predictor variables,
% builds observation matrices,
% and performs linear regression through least squares model identification

%%% INPUT
% data: original time series (N observations x M processes)
% jv: indexes of predicted series
% iv: indexes of predictors
% iv_lags: lags of predictors, in a row vector
% (they are all the same for the various predictors)

function out=lrp_idVAR(data,jv,iv,iv_lags)

[N,M]=size(data);

% this regression does not have constant term - remove the mean
% for m=1:M
%     data(:,m) = data(:,m) - mean(data(:,m));
% end

% compute the maxlag
maxlag=max(iv_lags);

My=data((maxlag+1:N)',jv); % observation matrix of the predicted variables
MX=[]; % observation matrix of the predictors
for l=1:length(iv_lags)
    for k=1:length(iv)
        MX=[MX data(maxlag+1-iv_lags(l):N-iv_lags(l),iv(k))];
    end
end

%eA_tmp=inv(MX'*MX)*MX'*My; % coefficients (least squares)
eA=(MX'*MX)\(MX'*My);
%p=size(eA,1)/M;

eu=My-MX*eA; % residuals
es2u=cov(eu); % covariance of residuals
es2y=cov(My); % covariance of predicted variables
erho2xy=1-det(es2u)/det(es2y); % squared correlation

out.eA=eA;
out.eu=eu;
out.es2u=es2u;
out.es2y=es2y;
out.erho2=erho2xy;
out.My=My;
out.MX=MX;

end
