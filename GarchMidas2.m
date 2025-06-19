function [estParams,EstParamCov,Variance,LongRunVar,ShortRunVar,logL,Table,logLikeSum] = GarchMidas2(y,varargin)
%GarchMidas: Maximum likelihood estimation of GARCH-MIDAS
% Parse inputs and set defaults.
periodDefault = 22;
nlagDefault = 10;
callerName = 'GarchMidas';
parseObj = inputParser;
addParameter(parseObj,'X',[],@(x)validateattributes(x,{'numeric'},{'2d'},callerName));
addParameter(parseObj,'NumLags',nlagDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'Period',periodDefault,@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'EstSample',[],@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},callerName));
addParameter(parseObj,'RollWindow',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'LogTau',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Beta2Para',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Options',[],@(x)validateattributes(x,{'struct','optim.options.Fmincon'},{},callerName));
addParameter(parseObj,'Mu0',[],@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'Alpha0',0.05,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'Beta0',0.9,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'Theta0',0.1,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'Theta1',0.1,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'W0',5,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'W1',5,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'M0',0.01,@(x)validateattributes(x,{'numeric'},{'scalar'},callerName));
addParameter(parseObj,'Gradient',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'AdjustLag',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'ThetaM',false,@(x)validateattributes(x,{'numeric','logical'},{'binary','nonempty'},callerName));
addParameter(parseObj,'Params',[],@(x)validateattributes(x,{'numeric'},{'column'},callerName));
addParameter(parseObj,'ZeroLogL',[],@(x)validateattributes(x,{'numeric'},{'2d','integer'},callerName));
parse(parseObj,varargin{:});
Regressor = parseObj.Results.X;   % X
period = parseObj.Results.Period;
nlag = parseObj.Results.NumLags;
estSample = parseObj.Results.EstSample;
rollWindow = parseObj.Results.RollWindow;
logTau = parseObj.Results.LogTau;
beta2Para = parseObj.Results.Beta2Para;
options = parseObj.Results.Options;
mu0 = parseObj.Results.Mu0;
alpha0 = parseObj.Results.Alpha0;
beta0 = parseObj.Results.Beta0;
theta0 = parseObj.Results.Theta0;
theta1 = parseObj.Results.Theta1;
w0 = parseObj.Results.W0;
w1 = parseObj.Results.W1;
m0 = parseObj.Results.M0;
autodiffFlag = parseObj.Results.Gradient;
adjustLag = parseObj.Results.AdjustLag;
thetaM = parseObj.Results.ThetaM;
estParams = parseObj.Results.Params;
EstParamCov = [];
zeroLogL = parseObj.Results.ZeroLogL;

if ~isempty(Regressor) %  如果Regressor不是空集
    if numel(Regressor(:,1)) ~= numel(y)  %  n = numel(A)  返回数组A中元素个数。
        error('Macroeconomic regressor must be a vector of the same length as y.')
    end
    
    if ~thetaM
        warning('Reset the Name-value pair ''thetaM'' to true, so that theta and M could be negative.')
        thetaM = true;
    end
end

% Replace missing values by the sample average 缺失值用样本平均值代替
y(isnan(y)) = nanmean(y);

% Divide data as the estimation and forecast samples将数据分为估计和预测样本
if isempty(estSample)  %如果estSample为空，则为全样本估计。
    estSample = numel(y); 
end
yFull = y;
y = y(1:estSample);
nobs = numel(y);  % nobs收益率样本内的个数
nMonth = ceil(nobs/period);  % nobs样本内的个数/period（22）=月份数

if ~isempty(Regressor)  % 如果Regressor不为空
%     Regressor(isnan(Regressor)) = nanmean(Regressor);
    RegressorFull = Regressor;
    Regressor = Regressor(1:estSample,:);
end

% Reshape the observation vector as a matrix, in which each column contains 
% observations in a week/month/quarter/year
Y = NaN(period,nMonth);  % 定义一个22行，月份数为列的NaN矩阵
Y(1:nobs) = y(:);     % nobs收益率样本内的个数

% Compute realized volatility or macroeconomic regressors
% Refer to Eq (6) in Engle et al. (2013)
if rollWindow
    
    % Realized volatility varies every period due to rolling windows 由于滚动窗口，每段时期的已实现波动都变化
    RV = zeros(nobs,2);   % 定义RV为一个样本内长度的序列
    for t = period:nobs   % 22：样本内长度
        if isempty(Regressor) % 如果Regressor是空集，那么RV用y^2代替
            RV(t,1) = nansum(y(t-period+1:t).^2); 
            RV=RV(t,1);
        else
%             RV(t,1) = nanmean(Regressor(t-period+1:t,1)); 
            RV(t,1) = nansum(Regressor(t-period+1:t,1).^2); 
            RV(t,2) = nanmean(Regressor(t-period+1:t,2)); 
        end
    end
    RV(1:period-1,1) = RV(period,1);
    RV(1:period-1,2) = RV(period,2);
else % 不滚动
    % Realized volatility is fixed in a week/month/quarter/year
    if isempty(Regressor)
        RV = nansum(Y.^2,1);
        % Odd days of the last week/month/quarter/year
        if mod(nobs,period) > 0  % mod(x,y):求模
            RV(end) = nansum(Y(end-period+1:end).^2);
        end
    else
        RegMat1 = NaN(period,nMonth);
        RegMat2 = NaN(period,nMonth);
        RegMat1(1:nobs) = nansum(Y.^2,1);
        RegMat2(1:nobs) = Regressor(:,2);
        RV(:,1) = nanmean(RegMat1,1);
        RV(:,2) = nanmean(RegMat2,1);
      
        % Odd days of the last week/month/quarter/year
        if mod(nobs,period) > 0
            RV(end,1) = nansum(Y(end-period+1:end).^2);
            RV(end,2) = nanmean(Regressor(end-period+1:end,2));
        end
    end
end

% Parameter estimation by maximum likelihood
if isempty(estParams)
    
    if isempty(mu0)
        mu0 = mean(y);
    end
    
    % Bounds for numerical optimization
    if beta2Para
        params0 = [mu0;alpha0;beta0;theta0;theta1;w0;w0;w1;w1;m0];
        if thetaM
            lb = [-Inf 0 0 -Inf -Inf 1.001  1.001 1.001  1.001  -Inf];
            ub = [ Inf 1 1  Inf Inf   50    50      50     50    Inf];
        else
            lb = [-Inf 0 0   0   0  1.001  0     1.001   0   0];
            ub = [ Inf 1 1  Inf Inf  50    50      50    50 Inf];
        end
    else
        params0 = [mu0;alpha0;beta0;theta0;theta1;w0;w1;m0];
        if thetaM
            lb = [-Inf 0 0 -Inf -Inf 1.001  1.001  -Inf];
            ub = [ Inf 1 1  Inf Inf    50    50     Inf];
        else
            lb = [-Inf 0 0   0    0  1.001  1.001   0];
            ub = [ Inf 1 1  Inf  Inf  50     50    Inf];
        end
        
    end
    
    % Objective function for maximum likelihood estimation
    if autodiffFlag
        if isempty(options)
            % options = optimoptions('fmincon','GradObj','on','DerivativeCheck','on','FinDiffType','central','Algorithm','interior-point','Display','notify-detailed');
            options = optimoptions('fmincon','GradObj','on','Algorithm','interior-point','Display','notify-detailed');
        end
        myfun = @(params)fNegMLGrad(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow);        
    else
        if isempty(options)
            options = optimoptions('fmincon','Algorithm','interior-point','Display','notify-detailed');
        end
        myfun = @(params) - nansum(fML(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow));        
    end
    
    % Use FMINCON for numerical optimization
    try
        estParams = fmincon(myfun,params0,[],[],[],[],lb,ub,[],options);
    catch Exception
        warning('FMINCON failed... Switch to FMINSEARCH. Error Message: %s',Exception.message)
        estParams = fminsearch(myfun,params0);
    end
    
    % Compute MLE covariance matrix
    [logL,Gradient] = fMLGrad(estParams,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow);
    BHHH = Gradient' * Gradient;
    if rcond(BHHH) < 1e-12
        warning('Covariance matrix of estimators cannot be computed precisely due to inversion difficulty.')        
    end
    EstParamCov = inv(BHHH);
    
    % Compute standard errors
    se2 = diag(EstParamCov);
    if any(se2<0)
        se2(se2<0) = NaN;
    end
    se = sqrt(se2);
    zstat = estParams(:) ./ se;
    pval = 0.5 * erfc(0.7071 * abs(zstat)) * 2;
    pval(pval<1e-6) = 0;
    
    % AIC and BIC
    logLikeSum = nansum(logL);
    numParam = length(params0);
    aic = -2*logLikeSum + 2*numParam;
    bic = -2*logLikeSum + numParam.*log(nobs);
    adjustSampleSize = sum(logL~=0 & ~isnan(logL));
    
    % Display the estimation results
    if autodiffFlag
        fprintf('Method: Maximum likelihood (Gradient)\n');
    else
        fprintf('Method: Maximum likelihood\n');
    end
    fprintf('Sample size: %d\n',nobs);
    fprintf('Adjusted sample size: %d\n',adjustSampleSize);
    fprintf('Logarithmic  likelihood: %12.6g\n',logLikeSum);
    fprintf('Akaike   info criterion: %12.6g\n',aic);
    fprintf('Bayesian info criterion: %12.6g\n',bic);
    columnNames = {'Coeff','StdErr','tStat','Prob'};
    if beta2Para
        rowNames = {'mu';'alpha';'beta';'theta1';'theta2';'w1-1';'w1-2';'w2-1';'w2-2';'m'};
    else
        rowNames = {'mu';'alpha';'beta';'theta1';'theta2';'w1';'w2';'m'};
    end
    try
        Table = table(estParams,se,zstat,pval,'RowNames',rowNames,'VariableNames',columnNames);
    catch
        Table = [estParams,se,zstat,pval];
    end
    disp(Table)
end

% Conditional variance with short-run and long-run components
nobsBig = numel(yFull);
nMonthBig = ceil(nobsBig/period);
Ybig = NaN(period,nMonthBig);
Ybig(1:nobsBig) = yFull(:);
if rollWindow
    RVBig = zeros(nobsBig,2);
    for t = period:nobsBig
        if isempty(Regressor)
            RVBig(t,1) = nansum(yFull(t-period+1:t).^2); 
            RVBig(t)=RVBig(t,1)
        else
            RVBig(t,1) = nansum(yFull(t-period+1:t).^2); 
            RVBig(t,2) = nanmean(RegressorFull(t-period+1:t,2)); 
        end
    end
    RVBig(1:period-1) = RVBig(period);
else
    if isempty(Regressor)
        RVBig = nansum(Ybig.^2,1);
        if mod(nobsBig,period) > 0
            RVBig(end) = nansum(yFull(end-period+1:end).^2);
        end        
    else
        RegMatBig1 = NaN(period,nMonthBig);
        RegMatBig2 = NaN(period,nMonthBig);
        RegMatBig1(1:nobs) = nansum(Ybig.^2,1);
        RegMatBig2(1:nobs) = RegressorFull(:,2);
        RVBig1 = nanmean(RegMatBig1,1);
        RVBig2 = nanmean(RegMatBig2,1);
        RVBig=[RVBig1,RVBig2];
        if mod(nobsBig,period) > 0
%             RVBig(end,1) = nanmean(RegressorFull(end-period+1:end,1));
            RVBig(end,1) = nansum(yFull(end-period+1:end).^2);
            RVBig(end,2) = nanmean(RegressorFull(end-period+1:end,2));
        end
    end
end
[~,Variance,ShortRunVar,LongRunVar] = fML(estParams,Ybig,RVBig,nlag,nobsBig,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow);

% One-step-ahead in-sample forecast validation
RealizedY2 = (yFull - estParams(1)).^2;
forecastError = Variance - RealizedY2;
estSampleRMSE = sqrt(mean(forecastError(1:estSample).^2));
outSampleRMSE = sqrt(mean(forecastError(estSample+1:end).^2)); 
fprintf('RMSE of one-step variance forecast (period 1 to %d): %5.3e.\n',estSample,estSampleRMSE);
if estSample < nobsBig
    fprintf('RMSE of one-step variance forecast (period %d to %d): %5.3e.\n',estSample+1,nobsBig,outSampleRMSE);
end
disp(' ')

end


%-------------------------------------------------------------------------
% Subfunction: MIDAS beta polynomial weights
function weights = midasBetaWeights(nlag,param1,param2)
seq = nlag:-1:1; 
if isempty(param2)    
    weights = (1-seq./nlag+10*eps).^(param1-1);    
else
    weights = (1-seq./nlag+10*eps).^(param1-1) .* (seq./nlag).^(param2-1);    
end
weights = weights ./ nansum(weights);
end


%-------------------------------------------------------------------------
% Subfunction: MIDAS beta polynomial weights
function weights = midasBetaWeightsCapital(nlag,param1,param2)
seq = nlag:-1:1;
if isempty(param2)
    weights = POWER(1-seq./nlag+10*eps, param1-1);
else
    weights = TIMES( POWER(1-seq./nlag+10*eps, param1-1), POWER(seq./nlag, param2-1)) ;
end
weights = RDIVIDE(weights,nansum(weights));
end


%-------------------------------------------------------------------------
% Subfunction: Compute log likelihood of GARCH-MIDAS
%              Also return conditional variance as a by-product
function [logL,Variance,ShortRun,LongRun] = fML(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow)

% Input Arguments:
%
% params: parameter vector [mu;alpha;beta;theta;w;m] 
% Y:      freq-by-nMonth matrix form of observations (nobs = freq * nMonth)
% RV:     1-by-nMonth realized volatility or macroeconomic variables
% nlag:   number of MIDAS lags
% nobs:   number of observations
% adjustLag: logical value on MIDAS lag adjustment for initial observations
% logTau: logical value on logarithmic long-run volatility component
% thetaM: logical value on not taking squares of theta and M
% zeroLogL:  a vector of indices that force zeros of the likelihood values 
% beta2Para: logical value on two parameter MIDAS Beta polynomial
% rollWindow: Logical value on rolling window long-run component
%
% Output Arguments:
%
% logL       nobs-by-1 log likelihood
% Variance   nobs-by-1 conditional variance
% ShortRun   nobs-by-1 short-run GARCH(1,1) component of conditional variance
% LongRun    nobs-by-1 long-run MIDAS component of conditional variance
%            LongRun values are fixed within a period
%

% Matrix dimension
[period,nMonth] = size(Y);

% Allocate parameters
mu0 = params(1);
alpha0 = params(2);
beta0 = params(3);
theta0 = params(4);
theta1 = params(5);
w10 = params(6);
if beta2Para
    w20 = params(7);
    w12 = params(8);
    w22 = params(9);    
    m0 = params(10);
else
    w20 = [];
    w22 = [];
    w12 = params(7);
    m0 = params(8);
end

% GARCH positive constraint
intercept = 1 - alpha0 - beta0;
if intercept < 0 || alpha0 < 0 || beta0 < 0 || alpha0 > 1 || beta0 > 1
    logL = -Inf(nobs,1);
    Variance = NaN(nobs,1);
    ShortRun = NaN(nobs,1);
    LongRun = NaN(nobs,1);
    return
end

% theta and m are squared for compatibility with others' codes
if ~thetaM
    theta0 = theta0 .* theta0;
    theta1 = theta1 .* theta1;
    m0 = m0 .* m0;
end

% Deflate observations and take squared residuals
Ydeflate = Y - mu0;
ResidSq = Ydeflate .* Ydeflate;

% Preallocate shortRun and Variance as a matrix, which will be reshaped
% Initial short-run component has the unconditional mean of one
% Initial long-run component has the unconditional mean of sample average
% Initialization will not affect likelihood computation
ShortRun = ones(period,nMonth);
if logTau
    tauAvg = exp(m0 + theta0 .* nanmean(RV(:,1))+ theta1 .* nanmean(RV(:,2)));
else
    tauAvg = m0+theta0 .* nanmean(RV(:,1))+ theta1 .* nanmean(RV(:,2));
end
Variance = tauAvg .* ones(period,nMonth);

% Conditional variance recursion
if rollWindow
    
    %---------------
    % Rolling Window
    %---------------
    
    % Compute MIDAS weights
    nlagBig = period*nlag;
    weights1 = midasBetaWeights(nlagBig,w10,w20);
    weights2 = midasBetaWeights(nlagBig,w12,w22);
    if adjustLag
        loopStart = 2;
    else
        loopStart = nlagBig + 1;
    end
    for t = loopStart:nobs
        
        % Compute long-run component
        % Refer to Eq (5) in Engle et al. (2013)
        if adjustLag && (t <= nlagBig+1)
            % Reduced MIDAS lags for initial observations
            weights1 = midasBetaWeights(t-1,w10,w20);
            weights2 = midasBetaWeights(t-1,w12,w22);
            tau = m0 + theta0 .* (weights1 * RV(1:t-1,1))+ theta1 .* (weights2 * RV(1:t-1,2));
        else            
            tau = m0 + theta0 .* (weights1 * RV(t-nlagBig:t-1,1))+ theta1 .* (weights2 * RV(t-nlagBig:t-1,2));
        end
      
        if logTau
            tau = exp(tau);
        end
      
        alphaTau = alpha0 ./ tau;
        
        % Compute short-run component
        % Refer to Eq (4) in Engle et al. (2013)
        ShortRun(t) = intercept + alphaTau .* ResidSq(t-1) + beta0 .* ShortRun(t-1);
        
        % Compute conditional variance
        % Refer to Eq (3) in Engle et al. (2013)
        Variance(t) = tau .* ShortRun(t);
        
    end    
    
else
    
    %---------------------------------------
    % Fixed tau in a week/month/quarter/year
    %---------------------------------------
    
    if adjustLag
        
        % Compute GARCH-MIDAS long-run and short-run variance components
        % The first column is unassigned due to missing presample values
        
        for t = 2:nMonth
            % Compute MIDAS weights
            % Reduced MIDAS lags for initial observations
            if t <= nlag+1
                weights1 = midasBetaWeights(t-1,w10,w20)';
                weights2 = midasBetaWeights(t-1,w12,w22)';
                RVuse = RV(1:t-1,:);
            else
                RVuse = RV(t-nlag:t-1,:);
            end
            
            % Compute long-run component
            % Refer to Eq (5) in Engle et al. (2013)
            if logTau
                tau = exp(m0 + theta0 .* (RVuse(:,1) * weights1)+ theta1 .* (RVuse(:,2) * weights2));
            else
                tau = m0 + theta0 .* (RVuse(:,1) * weights1)+ theta1 .* (RVuse(:,2) * weights2);
            end
            alphaTau = alpha0 ./ tau;
            
            % Compute short-run component
            % Refer to Eq (4) in Engle et al. (2013)
            for n = 1:period
                ind = (t-1)*period + n;
                ShortRun(ind) = intercept + alphaTau .* ResidSq(ind-1) + beta0 .* ShortRun(ind-1);
            end
            
            % Compute conditional variance
            % Refer to Eq (3) in Engle et al. (2013)
            Variance(:,t) = tau .* ShortRun(:,t);
        end
        
    else
        
        % Compute MIDAS weights
        weights1 = midasBetaWeights(nlag,w10,w20)';
        weights2 = midasBetaWeights(nlag,w12,w22)';
        % Compute GARCH-MIDAS long-run and short-run variance components
        % The first nlag columns are unassigned due to missing presample values
        for t = nlag+1:nMonth
            
            % Compute long-run component
            % Refer to Eq (5) in Engle et al. (2013)
            RVuse = RV(t-nlag:t-1);
            if logTau
                tau = exp(m0 + theta0 .* (RVuse(:,1) * weights1)+ theta1 .* (RVuse(:,2) * weights2));
            else
                tau = m0 + theta0 .* (RVuse(:,1) * weights1)+ theta1 .* (RVuse(:,2) * weights2);
            end
            alphaTau = alpha0 ./ tau;
            
            % Compute short-run component
            % Refer to Eq (4) in Engle et al. (2013)
            for n = 1:period
                ind = (t-1)*period + n;
                ShortRun(ind) = intercept + alphaTau .* ResidSq(ind-1) + beta0 .* ShortRun(ind-1);
            end
            
            % Compute conditional variance
            % Refer to Eq (3) in Engle et al. (2013)
            Variance(:,t) = tau .* ShortRun(:,t);
        end
    end
end

if any(Variance(:)<0)
    logL = -Inf(nobs,1);
    Variance = NaN(nobs,1);
    ShortRun = NaN(nobs,1);
    LongRun = NaN(nobs,1);
    % warning('Conditional variance falls below zero. Consider a larger ''m0'' name-value pair.')
    return
end

% Compute GARCH-MIDAS log likelihood 
logLMatrix = -0.5 .* ( log(2*pi.*Variance) + ResidSq ./ Variance);

% Presample log likelihood adjustment
if adjustLag
    logLMatrix(:,1) = 0;
else
    logLMatrix(:,1:nlag) = 0;    
end
logL = reshape(logLMatrix(1:nobs),nobs,1);

% User specified zero logL
logL(zeroLogL) = 0;

% If logMatrix are all NaNs, logL would be all zeros
if all(logL == 0) || nansum(logL) == 0
    logL = -Inf(nobs,1);
end    

% Reshape conditional variance as a column vector
if nargout > 1
    Variance = reshape(Variance(1:nobs),nobs,1);
    ShortRun = reshape(ShortRun(1:nobs),nobs,1);
    LongRun = Variance ./ ShortRun;
end

end


%-------------------------------------------------------------------------
% Subfunction: Compute log likelihood (vector) of GARCH-MIDAS
function logL = fMLCapital(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow)

% Matrix dimension 
[period,nMonth] = size(Y); 
 
% Allocate parameters 
mu0 = params(1); 
alpha0 = params(2); 
beta0 = params(3); 
theta0 = params(4);
theta1 = params(5);
w10 = params(6); 
if beta2Para
    w20 = params(7);
    w12 = params(8);
    w22 = params(9);    
    m0 = params(10);
else
    w20 = [];
    w22 = [];
    w12 = params(7);
    m0 = params(8);
end  
% GARCH positive constraint
intercept = 1 - alpha0 - beta0;
if intercept < 0 || alpha0 < 0 || beta0 < 0 || alpha0 > 1 || beta0 > 1
    logL = -Inf(nobs,1);
    return 
end 
 
% theta and m are squared for compatibility with others' codes 
if ~thetaM 
    theta0 = TIMES(theta0,theta0); 
    theta1 = TIMES(theta1,theta1);
    m0 = TIMES(m0,m0); 
end 
 
% Deflate observations and take squared residuals 
Ydeflate = Y - mu0; 
ResidSq = TIMES(Ydeflate,Ydeflate); 
 
% Preallocate shortRun and Variance as a matrix, which will be reshaped 
% Initial short-run component has the unconditional mean of one 
% Initial long-run component has the unconditional mean of sample average 
% Initialization will not affect likelihood computation 
ShortRun = ones(period,nMonth);   %生成全1矩阵
if logTau 
    tauAvg = EXP(m0 + theta0 .* nanmean(RV(:,1))+ theta1 .* nanmean(RV(:,2))); 
else 
    tauAvg = m0 + theta0 .* nanmean(RV(:,1))+ theta1 .* nanmean(RV(:,2)); 
end 
Variance = tauAvg .* ones(period,nMonth); 
 
% Conditional variance recursion
if rollWindow
    
    %---------------
    % Rolling Window
    %---------------
    
    % Compute MIDAS weights
    nlagBig = period*nlag;
    weights1 = midasBetaWeightsCapital(nlagBig,w10,w20);
    weights2 = midasBetaWeightsCapital(nlagBig,w12,w22);
    if adjustLag
        loopStart = 2;
    else
        loopStart = nlagBig + 1;
    end
    for t = loopStart:nobs
        
        % Compute long-run component
        % Refer to Eq (5) in Engle et al. (2013)
        if adjustLag && (t <= nlagBig+1)
            % Reduced MIDAS lags for initial observations
            weights1 = midasBetaWeightsCapital(t-1,w10,w20);
            weights2 = midasBetaWeightsCapital(t-1,w12,w22);
            tau = m0 + TIMES(theta0,MTIMES(weights1,RV(1:t-1,1)))+ TIMES(theta1,MTIMES(weights2,RV(1:t-1,2))); 
                 
        else            
            tau = m0 + TIMES(theta0,MTIMES(weights1,RV(t-nlagBig:t-1,1)))+ TIMES(theta1,MTIMES(weights2,RV(t-nlagBig:t-1,2))); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        if logTau
            tau = EXP(tau);
        end
        alphaTau = RDIVIDE(alpha0,tau);
        
        % Compute short-run component
        % Refer to Eq (4) in Engle et al. (2013)
        ShortRun(t) = intercept + TIMES(alphaTau,ResidSq(t-1)) + TIMES(beta0,ShortRun(t-1));
                
        % Compute conditional variance
        % Refer to Eq (3) in Engle et al. (2013)
        Variance(t) = TIMES(tau, ShortRun(t));
        
    end
    
else
    
    %---------------------------------------
    % Fixed tau in a week/month/quarter/year
    %---------------------------------------
    if adjustLag
        
        % Compute GARCH-MIDAS long-run and short-run variance components
        % The first column is unassigned due to missing presample values
        
        for t = 2:nMonth
            % Compute MIDAS weights
            % Reduced MIDAS lags for initial observations
            if t <= nlag+1
                weights1 = midasBetaWeightsCapital(t-1,w10,w20).';
                weights2 = midasBetaWeightsCapital(t-1,w12,w22).';
                RVuse = RV(1:t-1,:);
            else
                RVuse = RV(t-nlag:t-1,:);
            end
            
            % Compute long-run component
            % Refer to Eq (5) in Engle et al. (2013)
            if logTau
                tau = EXP(m0 + TIMES(theta0,MTIMES(RVuse(:,1),weights1))+ TIMES(theta1,MTIMES(RVuse(:,2),weights2)));
            else
                tau = m0 + TIMES(theta0,MTIMES(RVuse(:,1),weights1))+ TIMES(theta1,MTIMES(RVuse(:,2),weights2));
            end
            alphaTau = RDIVIDE(alpha0,tau);
            
            % Compute short-run component
            % Refer to Eq (4) in Engle et al. (2013)
            for n = 1:period
                ind = (t-1)*period + n;
                ShortRun(ind) = intercept + TIMES(alphaTau,ResidSq(ind-1)) + TIMES(beta0,ShortRun(ind-1));
            end
            
            % Compute conditional variance
            % Refer to Eq (3) in Engle et al. (2013)
            Variance(:,t) = TIMES(tau,ShortRun(:,t));
        end
        
    else
        
        % Compute MIDAS weights
        weights1 = midasBetaWeightsCapital(nlag,w10,w20).';
        weights2 = midasBetaWeightsCapital(nlag,w12,w22).';
        % Compute GARCH-MIDAS long-run and short-run variance components
        % The first nlag columns are unassigned due to missing presample values
        for t = nlag+1:nMonth
            
            % Compute long-run component
            % Refer to Eq (5) in Engle et al. (2013)
            RVuse = RV(t-nlag:t-1,:);
            if logTau
                tau = EXP(m0 + TIMES(theta0,MTIMES(RVuse(:,1),weights1))+ TIMES(theta1,MTIMES(RVuse(:,2),weights2)));
            else
                tau = m0 + TIMES(theta0,MTIMES(RVuse(:,1),weights1))+ TIMES(theta1,MTIMES(RVuse(:,2),weights2));
            end
            alphaTau = RDIVIDE(alpha0,tau);
            
            % Compute short-run component
            % Refer to Eq (4) in Engle et al. (2013)
            for n = 1:period
                ind = (t-1)*period + n;
                ShortRun(ind) = intercept + TIMES(alphaTau,ResidSq(ind-1)) + TIMES(beta0,ShortRun(ind-1));
            end
            
            % Compute conditional variance
            % Refer to Eq (3) in Engle et al. (2013)
            Variance(:,t) = TIMES(tau,ShortRun(:,t));
        end
    end
end
 
if any(Variance(:)<0) 
    logL = -Inf(nobs,1);      
    return 
end 
 
% Compute GARCH-MIDAS log likelihood 
logLMatrix = -0.5 .* (LOG(2*pi.*Variance) + RDIVIDE(ResidSq,Variance)); 
 
% Presample log likelihood adjustment 
if adjustLag 
    logLMatrix(:,1) = 0; 
else 
    logLMatrix(:,1:nlag) = 0; 
end 
logL = reshape(logLMatrix(1:nobs),nobs,1); 

% User specified zero logL
logL(zeroLogL) = 0;
 
% If logMatrix are all NaNs, logL would be all zeros 
if all(logL == 0) || nansum(logL) == 0
    logL = -Inf(nobs,1); 
end

end 


%-------------------------------------------------------------------------
% Subfunction: Compute log likelihood and gradient
function [logL,Gradient] = fMLGrad(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow)

logL = fML(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow);
if nargout == 1    
    return
end

fun = @(params) fMLCapital(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow);
x = params;

nfunVal = numel(logL);
nparams = numel(x);
fval = zeros(nfunVal,nparams);
Gradient = zeros(nfunVal,nparams);

for m = 1:nparams
    
    % One element of x carries an imaginary number
    xComplex = x;
    xComplex(m) = xComplex(m) + 1i;
    
    % The real component is the function value 
    % The imaginary component is the derivative
    complexValue = fun(xComplex);
    fval(:,m) = real(complexValue);
    Gradient(:,m) = imag(complexValue);
end

if any(norm(bsxfun(@minus,fval, logL))>1e-6)
    warning('Automatic gradients might be falsely computed. Capitalize arithmetic operator and utility function names.')
end

end


%-------------------------------------------------------------------------
% Subfunction: Compute negative log likelihood (scalar) and gradient
function [logL,Gradient] = fNegMLGrad(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow)
[logL,Gradient] = fMLGrad(params,Y,RV,nlag,nobs,adjustLag,logTau,thetaM,zeroLogL,beta2Para,rollWindow);
logL = - nansum(logL,1);
Gradient = - nansum(Gradient,1);
end


%-------------------------------------------------------------------------
function B = EXP(A)
Areal = real(A);
Aimag = imag(A);
Breal = exp(Areal);
Bimag = Breal .* Aimag;
B = complex(Breal, Bimag);
end


%-------------------------------------------------------------------------
function B = LOG(A)
Areal = real(A);
Aimag = imag(A);
Breal = log(Areal);
Bimag = Aimag ./ Areal;
B = complex(Breal, Bimag);
end


%-------------------------------------------------------------------------
function C = MTIMES(A,B)
C = A * B + imag(A) * imag(B);
end


%-------------------------------------------------------------------------
function C = POWER(A,B)
if isscalar(B) && (B == 3)
    C = TIMES(TIMES(A, A), A);
elseif isscalar(B) && (B == 2)
    C = TIMES(A, A);
elseif isscalar(B) && (B == 1)
    C = A;
elseif isscalar(B) && (B == 0.5)
    C = SQRT(A);
elseif isscalar(B) && (B == 0)
    C = complex(ones(size(A)),0);
elseif isscalar(B) && (B == -0.5)
    C = RDIVIDE(1, SQRT(A));
elseif isscalar(B) && (B == -1)
    C = RDIVIDE(1, A);
elseif isscalar(B) && (B == -2)
    C = RDIVIDE(1, TIMES(A,A));
elseif all(A(:)>0)
    C = EXP(TIMES(B, LOG(A)));
else
    % This is not necessarily correct
    Creal = real(A).^real(B);
    Cimag = imag(exp(B.*log(abs(A))));
    C = complex(Creal,Cimag);
end
end


%-------------------------------------------------------------------------
function C = RDIVIDE(A,B)
Areal = real(A);  %real是求复数的实数部分函数
Aimag = imag(A);
Breal = real(B);
Bimag = imag(B);
Creal = Areal ./ Breal;
Cimag = (Aimag .* Breal - Areal .* Bimag) ./ (Breal.*Breal);
C = complex(Creal, Cimag);
end


%-------------------------------------------------------------------------
function C = TIMES(A,B)
C = A .* B + imag(A) .* imag(B);   %imag是求复数的虚数部分函数
end

