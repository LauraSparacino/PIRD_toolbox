%% STATIC LINEAR REGRESSION MODEL
%% given data matrix and indexes of predicted and predictor variables, builds observation matrices, and performs linear regression through least squares model identification

%%% INPUT
% data: original time series
% jv: indexes of predicted series
% iv: indexes of predictors
% iv_lags: lagsof predictors, in a row vector (they are all the same for the various predictors)

function out=lrp_LinRegStatic(data,jv,iv)

[N,Q]=size(data);

% this regression does not have constant term - remove the mean
for m=1:Q
    data(:,m) = data(:,m) - mean(data(:,m));
end

MY=data(:,jv); % observation matrix of the predicted variables
MX=data(:,iv); % observation matrix of the predictor variables

% eA2=inv(MX'*MX)*MX'*MY; % coefficients (least squares)
eA=(MX'*MX)\MX'*MY; % coefficients (least squares)
eu=MY-MX*eA; % residuals
es2u=cov(eu); %covariance of residuals
es2y=cov(MY); %covariance of predicted variables
erho2xy=1-det(es2u)/det(es2y); %squared correlation

% out.eA2=eA2;
out.eA=eA;
out.eu=eu;
out.es2u=es2u;
out.es2y=es2y;
out.erho2=erho2xy;
% out.My=MY;
% out.MX=MX;

end
