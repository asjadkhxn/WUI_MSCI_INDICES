%%  ------Table 2. In-sample estimation results of GARCH-MIDAS models.----------
clc;
clear;

% =============== 1) Load Data ===============
RTD = xlsread('RTD_GPR_V.xlsx','RTD','C2:C4230');
nobs = size(RTD,1);

LogTau = 1;
thetaM = 1;
Beta2Para = 0; 
outsample = 0;  
Rollwindow = 1;  

xMonth = xlsread('RTD_GPR_V.xlsx','GPR','C3:C205')./100 ; 

% Repeat the monthly value throughout the days in that month
[~, yDate] = xlsread('RTD_GPR_V.xlsx', 'RTD','A2:A4230');
[~, yDateMonth] = datevec(yDate);
xDay = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)    
        count = count + 1;
        if count > length(xMonth)
            break
        end
    end
    xDay(t) = xMonth(count);
end

period = 22;
numLags = 36;

% =============== 2) RV Model ===============
[estParams_RV, EstParamCov_RV, Variance_RV, LongRunVar_RV] = GarchMidas(RTD, ...
    'Period',period,'NumLags',numLags,'ThetaM',1,'RollWindow',1,'estSample',4229);

y2 = NaN(nobs,1);
for t = 1:nobs
    y2(t) = estParams_RV(1) + sqrt(Variance_RV(t));
end
residual_RV = RTD - y2;

% =============== 3) GPR Model ===============
[estParams_GPR, EstParamCov_GPR, Variance_GPR, LongRunVar_GPR] = GarchMidas(RTD, ...
    'Period',period,'NumLags',numLags,'X',xDay,'ThetaM',1,'RollWindow',1,'estSample',4229);

y1 = NaN(nobs,1);
for t = 1:nobs
    y1(t) = estParams_GPR(1) + sqrt(Variance_GPR(t));
end
residual_x = RTD - y1;

seq = (period*numLags+1:nobs)';
year = linspace(2007.5,2024.3,nobs);
condition_variance2 = sqrt(252*Variance_GPR(seq));
long_var2          = sqrt(252*LongRunVar_GPR(seq));
GPR_var            = [condition_variance2 long_var2];

% =============== 4) RV+GPR Model ===============
x21 = [RTD xDay];
[estParams_RVGPR, EstParamCov_RVGPR, Variance_RVGPR, LongRunVar_RVGPR] = GarchMidas2(RTD, ...
    'X',x21,'Period',period,'NumLags',numLags, ...
    'LogTau',LogTau, ...
    'EstSample',4229,'Rollwindow',Rollwindow,'thetaM',thetaM, ...
    'Beta2Para',Beta2Para);

y3 = NaN(nobs,1);
for t = 1:nobs
    y3(t) = estParams_RVGPR(1) + sqrt(Variance_RVGPR(t));
end
residual_x_RV_GPR = RTD - y3;

% =============== 5) Save Residuals (if needed) ===============
save('Table2_Residuals.mat', 'residual_RV', 'residual_x', 'residual_x_RV_GPR');

% =============== 6) Extract Parameter Estimates & Std. Errors ===============
% Rename variables to avoid conflicts:
RV_params    = estParams_RV(:);
RV_stdErr    = sqrt(diag(EstParamCov_RV));

GPR_params   = estParams_GPR(:);
GPR_stdErr   = sqrt(diag(EstParamCov_GPR));

RVGPR_params = estParams_RVGPR(:);
RVGPR_stdErr = sqrt(diag(EstParamCov_RVGPR));

% =============== 7) Build "Table 2" for Excel ===============
% Determine number of parameters based on the RV model:
numParams = length(RV_params);
% Default names for up to 8 parameters
defaultParamNames = {'\mu','\alpha','\beta','\theta_1','\theta_2','\omega_1','\omega_2','m'};
% Use only as many names as there are estimated parameters:
paramNames = defaultParamNames(1:numParams);

% Functions to format numbers:
formatParams = @(p) sprintf('%.4f', p);
formatStdErr = @(s) sprintf('(%.4f)', s);

% Create cell arrays for each model's formatted values:
rvEst  = cell(1, numParams);
rvSE   = cell(1, numParams);
gprEst = cell(1, numParams);
gprSE  = cell(1, numParams);
rvgprEst  = cell(1, numParams);
rvgprSE   = cell(1, numParams);

for i = 1:numParams
    rvEst{i}    = formatParams(RV_params(i));
    rvSE{i}     = formatStdErr(RV_stdErr(i));
    
    gprEst{i}   = formatParams(GPR_params(i));
    gprSE{i}    = formatStdErr(GPR_stdErr(i));
    
    rvgprEst{i} = formatParams(RVGPR_params(i));
    rvgprSE{i}  = formatStdErr(RVGPR_stdErr(i));
end

% Build the rows for the table. Each row has (numParams+1) cells.
headerRow   = [{' '}, paramNames];
rvRowEst    = [{'RV'}, rvEst];
rvRowSE     = [{''}, rvSE];
gprRowEst   = [{'GPR'}, gprEst];
gprRowSE    = [{''}, gprSE];
rvgprRowEst = [{'RV+GPR'}, rvgprEst];
rvgprRowSE  = [{''}, rvgprSE];

% Build a notes row with (numParams+1) cells:
notesRow = cell(1, numParams+1);
notesRow{1} = 'Notes:';
if numParams+1 >= 2
    notesRow{2} = 'Standard errors in parentheses, *** 1% significance.';
end
for i = 3:(numParams+1)
    notesRow{i} = '';
end

% Combine all rows into a single cell array:
table2Data = [headerRow;
              rvRowEst;
              rvRowSE;
              gprRowEst;
              gprRowSE;
              rvgprRowEst;
              rvgprRowSE;
              notesRow];

% Convert cell array to table and write to Excel without auto headers:
T = cell2table(table2Data);
writetable(T, 'Table2_InSampleResults_v5.xlsx', 'WriteVariableNames', false);


%%            Table 3: Out-of-sample forecasting results script

clc;
clear;

% 1) Load the data
RTD = xlsread('RTD_GPR_V.xlsx','RTD','C2:C4230');
nobs = size(RTD,1);

LogTau    = 1;
thetaM    = 1;
Beta2Para = 0;
Rollwindow = 1;

xMonth = xlsread('RTD_GPR_V.xlsx','GPR','C3:C205') ./ 100;

% Repeat the monthly value throughout the days in that month
[~, yDate] = xlsread('RTD_GPR_V.xlsx','RTD','A2:A4230');
[~, yDateMonth] = datevec(yDate);
xDay = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)
        count = count + 1;
        if count > length(xMonth)
            break;
        end
    end
    xDay(t) = xMonth(count);
end

% Common parameters
period  = 22;
numLags = 36;

% 2) Define the forecast horizons to evaluate
horizons = [300, 600, 1000, 1200];

% 3) Prepare a cell array to store all rows for Table 3
table3Data = {};

% 4) Loop over each horizon h
for h = horizons
    
    % ========== 4.1) Sub-sample size ==========
    outSample = h;
    EstSample = nobs - outSample;  % in-sample length
    
    % -------------------------------------------------------
    %    (A) RV Model
    % -------------------------------------------------------
    [estParams_RV, EstParamCov_RV, Variance_RV, LongRunVar_RV] = ...
        GarchMidas(RTD, ...
                   'Period',     period, ...
                   'NumLags',    numLags, ...
                   'ThetaM',     1, ...
                   'RollWindow', 1, ...
                   'estSample',  EstSample);

    % Forecast variance for out-of-sample period
    Forecast_variance_RV = Variance_RV(EstSample+1 : end);
    actual_RV            = RTD(EstSample+1 : end).^2;
    
    % Compute forecast metrics
    RMSE_RV = sqrt(mean((Forecast_variance_RV - actual_RV).^2, 'omitnan'));
    RMAE_RV = sqrt(mean(abs(Forecast_variance_RV - actual_RV), 'omitnan'));
    MSPE_RV = mean((sqrt(Forecast_variance_RV+0.0001) - sqrt(actual_RV+0.0001)).^2, 'omitnan');
    MAPE_RV = mean(abs(sqrt(Forecast_variance_RV+0.0001) - sqrt(actual_RV+0.0001)), 'omitnan');
    RMSD_RV = sqrt(mean((sqrt(Forecast_variance_RV+0.0001) - sqrt(actual_RV+0.0001)).^2, 'omitnan'));
    RMAD_RV = sqrt(abs(mean((sqrt(Forecast_variance_RV+0.0001) - sqrt(actual_RV+0.0001)).^2, 'omitnan')));

    % -------------------------------------------------------
    %    (B) GPR Model
    % -------------------------------------------------------
    [estParams_GPR, EstParamCov_GPR, Variance_GPR, LongRunVar_GPR] = ...
        GarchMidas(RTD, ...
                   'Period',     period, ...
                   'NumLags',    numLags, ...
                   'X',          xDay, ...
                   'ThetaM',     1, ...
                   'RollWindow', 1, ...
                   'estSample',  EstSample);

    Forecast_variance_GPR = Variance_GPR(EstSample+1 : end);
    actual_GPR            = RTD(EstSample+1 : end).^2;
    
    RMSE_GPR = sqrt(mean((Forecast_variance_GPR - actual_GPR).^2, 'omitnan'));
    RMAE_GPR = sqrt(mean(abs(Forecast_variance_GPR - actual_GPR), 'omitnan'));
    MSPE_GPR = mean((sqrt(Forecast_variance_GPR+0.0001) - sqrt(actual_GPR+0.0001)).^2, 'omitnan');
    MAPE_GPR = mean(abs(sqrt(Forecast_variance_GPR+0.0001) - sqrt(actual_GPR+0.0001)), 'omitnan');
    RMSD_GPR = sqrt(mean((sqrt(Forecast_variance_GPR+0.0001) - sqrt(actual_GPR+0.0001)).^2, 'omitnan'));
    RMAD_GPR = sqrt(abs(mean((sqrt(Forecast_variance_GPR+0.0001) - sqrt(actual_GPR+0.0001)).^2, 'omitnan')));

    % -------------------------------------------------------
    %    (C) RV + GPR Model
    % -------------------------------------------------------
    x21 = [RTD, xDay];
    [estParams_RVGPR, EstParamCov_RVGPR, Variance_RVGPR, LongRunVar_RVGPR] = ...
        GarchMidas2(RTD, ...
                    'X',         x21, ...
                    'Period',    period, ...
                    'NumLags',   numLags, ...
                    'LogTau',    LogTau, ...
                    'EstSample', EstSample, ...
                    'Rollwindow',Rollwindow, ...
                    'thetaM',    thetaM, ...
                    'Beta2Para',Beta2Para);

    Forecast_variance_RVGPR = Variance_RVGPR(EstSample+1 : end);
    actual_RVGPR            = RTD(EstSample+1 : end).^2;
    
    RMSE_RVGPR = sqrt(mean((Forecast_variance_RVGPR - actual_RVGPR).^2, 'omitnan'));
    RMAE_RVGPR = sqrt(mean(abs(Forecast_variance_RVGPR - actual_RVGPR), 'omitnan'));
    MSPE_RVGPR = mean((sqrt(Forecast_variance_RVGPR+0.0001) - sqrt(actual_RVGPR+0.0001)).^2, 'omitnan');
    MAPE_RVGPR = mean(abs(sqrt(Forecast_variance_RVGPR+0.0001) - sqrt(actual_RVGPR+0.0001)), 'omitnan');
    RMSD_RVGPR = sqrt(mean((sqrt(Forecast_variance_RVGPR+0.0001) - sqrt(actual_RVGPR+0.0001)).^2, 'omitnan'));
    RMAD_RVGPR = sqrt(abs(mean((sqrt(Forecast_variance_RVGPR+0.0001) - sqrt(actual_RVGPR+0.0001)).^2, 'omitnan')));

    % 4.2) Build a block of rows for the current horizon h
    %   (a) Heading row: "h=XXX" + column headers
    headingRow = {
        sprintf('h=%d', h), ...
        'RMSE','RMAE','MSPE','MAPE','RMSD','RMAD'
    };

    %   (b) 1 row per model
    rowRV = {
        'RV', ...
        sprintf('%.4f', RMSE_RV), ...
        sprintf('%.4f', RMAE_RV), ...
        sprintf('%.4f', MSPE_RV), ...
        sprintf('%.4f', MAPE_RV), ...
        sprintf('%.4f', RMSD_RV), ...
        sprintf('%.4f', RMAD_RV)
    };

    rowGPR = {
        'GPR', ...
        sprintf('%.4f', RMSE_GPR), ...
        sprintf('%.4f', RMAE_GPR), ...
        sprintf('%.4f', MSPE_GPR), ...
        sprintf('%.4f', MAPE_GPR), ...
        sprintf('%.4f', RMSD_GPR), ...
        sprintf('%.4f', RMAD_GPR)
    };

    rowRVGPR = {
        'RV+GPR', ...
        sprintf('%.4f', RMSE_RVGPR), ...
        sprintf('%.4f', RMAE_RVGPR), ...
        sprintf('%.4f', MSPE_RVGPR), ...
        sprintf('%.4f', MAPE_RVGPR), ...
        sprintf('%.4f', RMSD_RVGPR), ...
        sprintf('%.4f', RMAD_RVGPR)
    };

    %   (c) Append these 4 new rows to our table3Data cell array
    table3Data = [
        table3Data;
        headingRow;
        rowRV;
        rowGPR;
        rowRVGPR
    ];
end

% 5) Add a final "Notes" row (span as many columns as needed)
notesRow = {
    'Notes:', ...
    '(1) h is the short-term out-of-sample prediction horizon of 300, 600, 1000, 1200 days.', ...
    '(2) RMSE = root mean square error; RMAE = root mean absolute error;', ...
    '    MSPE = mean square percentage error; MAPE = mean absolute percentage error;', ...
    '    RMSD = root mean square deviation; RMAD = relative mean absolute deviation.', ...
    'The smaller, the better for all the above indicators.'
};
% If needed, ensure the row has exactly 7 columns:
numCols = 7;
if length(notesRow) < numCols
    notesRow{numCols} = '';
end

table3Data = [table3Data; notesRow];

% 6) Convert the final cell array to a table and write to Excel
T = cell2table(table3Data);
writetable(T, 'Table3_OutOfSampleResults_v2.xlsx', 'WriteVariableNames', false);

disp('=== Done! See "Table3_OutOfSampleResults_v2.xlsx" for your final table. ===');



%% Table 4: MCS tests results for multiple horizons

clc;
clear;

% 1) Add path if needed (where "mcs" and "mcs_resort" are located)
addpath('test_code');

% 2) Load Data
RTD   = xlsread('RTD_GPR_V.xlsx','RTD','C2:C4230');
nobs  = size(RTD,1);

LogTau     = 1;
thetaM     = 1;
Beta2Para  = 0;
Rollwindow = 1;

xMonth = xlsread('RTD_GPR_V.xlsx','GPR','C3:C205') ./ 100;

% Repeat the monthly value throughout the days in that month
[~, yDate] = xlsread('RTD_GPR_V.xlsx', 'RTD','A2:A4230');
[~, yDateMonth] = datevec(yDate);
xDay = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)
        count = count + 1;
        if count > length(xMonth)
            break
        end
    end
    xDay(t) = xMonth(count);
end

period  = 22;
numLags = 36;

% 3) Define horizons
horizons = [300, 600, 1000, 1200];

% 4) Prepare a cell array to hold the entire "Table 4"
table4Data = {};

% 5) Loop over each horizon h, run the MCS tests, and append the results
for h = horizons
    
    % ------------------------ 5.1) Subset / Estimation Sample ------------------------
    outSample = h;
    EstSample = nobs - outSample;
    
    % --------------- (A) RV Model ---------------
    [estParams,EstParamCov,Variance,LongRunVar] = ...
        GarchMidas(RTD,'Period',period,'NumLags',numLags,'ThetaM',1, ...
                   'RollWindow',1,'estSample',EstSample);
    Forecast_variance_RV = Variance(EstSample+1:end);
    
    % --------------- (B) GPR Model ---------------
    [estParams,EstParamCov,Variance,LongRunVar] = ...
        GarchMidas(RTD,'Period',period,'NumLags',numLags,'X',xDay,'ThetaM',1, ...
                   'RollWindow',1,'estSample',EstSample);
    Forecast_variance_GPR = Variance(EstSample+1:end);
    
    % --------------- (C) RV+GPR Model ---------------
    x21 = [RTD, xDay];
    [estParams,EstParamCov,Variance,LongRunVar] = ...
        GarchMidas2(RTD,'X',x21,'Period',period,'NumLags',numLags, ...
                    'LogTau',LogTau,'EstSample',EstSample, ...
                    'Rollwindow',Rollwindow,'thetaM',thetaM, ...
                    'Beta2Para',Beta2Para);
    Forecast_variance_RVGPR = Variance(EstSample+1:end);
    
    % --------------- (D) Compute Losses for MCS ---------------
    actual = RTD(EstSample+1:end).^2;
    FC     = [Forecast_variance_RV, Forecast_variance_GPR, Forecast_variance_RVGPR];
    N_FC   = size(FC,2);  % should be 3 models

    % QLIKE
    LOSS_QLIKE = log(FC) + repmat(actual,1,N_FC)./FC;  
    % MSE
    LOSS_MSE = (FC - repmat(actual,1,N_FC)).^2;
    % MAE
    LOSS_MAE = abs(FC - repmat(actual,1,N_FC));
    % HMSE
    LOSS_HMSE = (1 - FC./repmat(actual,1,N_FC)).^2;
    % HMAE
    LOSS_HMAE = abs(1 - FC./repmat(actual,1,N_FC));
    
    % --------------- (E) Run MCS for each measure ---------------
    % QLIKE
    [incR,pvR,excR,incSQ,pvSQ,excSQ] = mcs(LOSS_QLIKE,0.0001,10000,25);
    p_mcs_QLIKE    = mcs_resort(incR,pvR,excR);
    p_mcs_QLIKE_SQ = mcs_resort(incSQ,pvSQ,excSQ);

    % MSE
    [incR,pvR,excR,incSQ,pvSQ,excSQ] = mcs(LOSS_MSE,0.0001,10000,25);
    p_mcs_MSE    = mcs_resort(incR,pvR,excR);
    p_mcs_MSE_SQ = mcs_resort(incSQ,pvSQ,excSQ);

    % MAE
    [incR,pvR,excR,incSQ,pvSQ,excSQ] = mcs(LOSS_MAE,0.0001,10000,25);
    p_mcs_MAE    = mcs_resort(incR,pvR,excR);
    p_mcs_MAE_SQ = mcs_resort(incSQ,pvSQ,excSQ);

    % HMSE
    [incR,pvR,excR,incSQ,pvSQ,excSQ] = mcs(LOSS_HMSE,0.0001,10000,25);
    p_mcs_HMSE    = mcs_resort(incR,pvR,excR);
    p_mcs_HMSE_SQ = mcs_resort(incSQ,pvSQ,excSQ);

    % HMAE
    [incR,pvR,excR,incSQ,pvSQ,excSQ] = mcs(LOSS_HMAE,0.0001,10000,25);
    p_mcs_HMAE    = mcs_resort(incR,pvR,excR);
    p_mcs_HMAE_SQ = mcs_resort(incSQ,pvSQ,excSQ);

    % --------------- (F) Combine all p-values: 3 rows x 10 columns ---------------
    % Columns in p_mcs_ALL3:
    %  1: QLIKE_TR,   2: QLIKE_TSQ,
    %  3: MSE_TR,     4: MSE_TSQ,
    %  5: MAE_TR,     6: MAE_TSQ,
    %  7: HMSE_TR,    8: HMSE_TSQ,
    %  9: HMAE_TR,    10: HMAE_TSQ
    p_mcs_ALL3 = [
        p_mcs_QLIKE,    p_mcs_QLIKE_SQ, ...
        p_mcs_MSE,      p_mcs_MSE_SQ, ...
        p_mcs_MAE,      p_mcs_MAE_SQ, ...
        p_mcs_HMSE,     p_mcs_HMSE_SQ, ...
        p_mcs_HMAE,     p_mcs_HMAE_SQ
    ];
    
    % ------------------------ 5.2) Build the sub-table for horizon h ------------------------
    % (a) First heading row: "h=XXX" + measure names
    % We want a total of 11 columns:
    %   - 1 for "h=300"
    %   - 2 each for QLIKE, MSE, MAE, HMSE, HMAE (so 2*5 = 10) => total 11
    headingRow1 = {
        sprintf('h=%d', h), ...
        'QLIKE','', 'MSE','', 'MAE','', 'HMSE','', 'HMAE',''
    };
    
    % (b) Second heading row: blank in col 1, then "TR","TSQ" pairs
    headingRow2 = {
        '', ...
        'TR','TSQ','TR','TSQ','TR','TSQ','TR','TSQ','TR','TSQ'
    };

    % (c) 3 model rows: RV, GPR, RV+GPR
    % Convert numeric p-values into string format (4 decimals).
    % If you want to bold numbers > 0.10, you'd have to add special text or handle it in Excel manually.
    rowRV = {
        'RV', ...
        sprintf('%.4f', p_mcs_ALL3(1,1)), sprintf('%.4f', p_mcs_ALL3(1,2)), ...
        sprintf('%.4f', p_mcs_ALL3(1,3)), sprintf('%.4f', p_mcs_ALL3(1,4)), ...
        sprintf('%.4f', p_mcs_ALL3(1,5)), sprintf('%.4f', p_mcs_ALL3(1,6)), ...
        sprintf('%.4f', p_mcs_ALL3(1,7)), sprintf('%.4f', p_mcs_ALL3(1,8)), ...
        sprintf('%.4f', p_mcs_ALL3(1,9)), sprintf('%.4f', p_mcs_ALL3(1,10))
    };

    rowGPR = {
        'GPR', ...
        sprintf('%.4f', p_mcs_ALL3(2,1)), sprintf('%.4f', p_mcs_ALL3(2,2)), ...
        sprintf('%.4f', p_mcs_ALL3(2,3)), sprintf('%.4f', p_mcs_ALL3(2,4)), ...
        sprintf('%.4f', p_mcs_ALL3(2,5)), sprintf('%.4f', p_mcs_ALL3(2,6)), ...
        sprintf('%.4f', p_mcs_ALL3(2,7)), sprintf('%.4f', p_mcs_ALL3(2,8)), ...
        sprintf('%.4f', p_mcs_ALL3(2,9)), sprintf('%.4f', p_mcs_ALL3(2,10))
    };

    rowRVGPR = {
        'RV+GPR', ...
        sprintf('%.4f', p_mcs_ALL3(3,1)), sprintf('%.4f', p_mcs_ALL3(3,2)), ...
        sprintf('%.4f', p_mcs_ALL3(3,3)), sprintf('%.4f', p_mcs_ALL3(3,4)), ...
        sprintf('%.4f', p_mcs_ALL3(3,5)), sprintf('%.4f', p_mcs_ALL3(3,6)), ...
        sprintf('%.4f', p_mcs_ALL3(3,7)), sprintf('%.4f', p_mcs_ALL3(3,8)), ...
        sprintf('%.4f', p_mcs_ALL3(3,9)), sprintf('%.4f', p_mcs_ALL3(3,10))
    };

    % (d) Append these 5 new rows to table4Data
    table4Data = [
        table4Data;
        headingRow1;
        headingRow2;
        rowRV;
        rowGPR;
        rowRVGPR
    ];
    
end

% 6) Add a final "Notes" row. We have 11 columns total.
notesRow = {
    'Notes:', ...
    'The values in the table are p-values, with numbers > 0.1000 in bold, indicating the model performed significantly better.', ...
    '0.1000 represents the highest accuracy of the model''s predictions.', ...
    '','', '','', '','', '',''
};

% Ensure we have exactly 11 cells in notesRow
if length(notesRow) < 11
    notesRow{11} = '';
end

table4Data = [table4Data; notesRow];

% 7) Convert cell array to table and write to Excel/CSV
T = cell2table(table4Data);
writetable(T, 'Table4_MCS_Results.xlsx', 'WriteVariableNames', false);

disp('=== Done! Check "Table4_MCS_Results.xlsx" for Table 4. ===');

%% Table 5: Estimation results of GARCH-MIDAS-GPR and GARCH-MIDAS-GPRH


clc;
clear;

% =============== 1) Load Data ===============
RTD = xlsread('RTD_GPR_V.xlsx','RTD','C2:C4230');
nobs = size(RTD,1);

LogTau     = 1;
thetaM     = 1;
Beta2Para  = 0; 
Rollwindow = 1;

% Read the monthly regressor series from two sheets:
xMonth_GPR   = xlsread('RTD_GPR_V.xlsx','GPR','C3:C205')./100;
xMonth_GPRH  = xlsread('RTD_GPR_V.xlsx','GPRH','C3:C205')./100;

% =============== 2) Build Daily Series for Each Model ===============
[~, yDate] = xlsread('RTD_GPR_V.xlsx','RTD','A2:A4230');
[~, yDateMonth] = datevec(yDate);

% Build daily series for GPR:
xDay_GPR = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)
        count = count + 1;
        if count > length(xMonth_GPR)
            break;
        end
    end
    xDay_GPR(t) = xMonth_GPR(count);
end

% Build daily series for GPRH:
xDay_GPRH = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)
        count = count + 1;
        if count > length(xMonth_GPRH)
            break;
        end
    end
    xDay_GPRH(t) = xMonth_GPRH(count);
end

% =============== 3) Estimate GPR Model ===============
period  = 22;
numLags = 36;
estSample = 4229;  % adjust as needed

[estParams_GPR, EstParamCov_GPR, ~, ~] = ...
    GarchMidas(RTD, 'Period', period, 'NumLags', numLags, ...
               'X', xDay_GPR, 'ThetaM', 1, 'RollWindow', 1, 'estSample', estSample);

GPR_params = estParams_GPR(:);
GPR_stdErr = sqrt(diag(EstParamCov_GPR));

% =============== 4) Estimate GPRH Model ===============
% Replace with your actual GPRH estimation if different.
[estParams_GPRH, EstParamCov_GPRH, ~, ~] = ...
    GarchMidas(RTD, 'Period', period, 'NumLags', numLags, ...
               'X', xDay_GPRH, 'ThetaM', 1, 'RollWindow', 1, 'estSample', estSample);

GPRH_params = estParams_GPRH(:);
GPRH_stdErr = sqrt(diag(EstParamCov_GPRH));

% =============== 5) Format Table 5 ===============
% Define parameter names (adjust if the number of parameters differs)
paramNames = {'\mu','\alpha','\beta','\theta','\omega','m'};

% Formatting functions
formatParams = @(p) sprintf('%.4f', p);
formatStdErr = @(s) sprintf('(%.4f)', s);

numParams = length(paramNames);

% Convert GPR estimates to cell arrays of strings
gprEst = cell(1, numParams);
gprSE  = cell(1, numParams);
for i = 1:numParams
    gprEst{i} = formatParams(GPR_params(i));
    gprSE{i}  = formatStdErr(GPR_stdErr(i));
end

% Convert GPRH estimates to cell arrays of strings
gprhEst = cell(1, numParams);
gprhSE  = cell(1, numParams);
for i = 1:numParams
    gprhEst{i} = formatParams(GPRH_params(i));
    gprhSE{i}  = formatStdErr(GPRH_stdErr(i));
end

% Build the 2D cell array. Each row must have exactly 1 + numParams (7) columns.
table5Data = {
    % Row 1: Header row
    ' ',       paramNames{:};          % 1 + 6 = 7 columns
    % Row 2: GPR estimates
    'GPR',     gprEst{:};              % 7 columns
    % Row 3: GPR standard errors
    '',        gprSE{:};               % 7 columns
    % Row 4: GPRH estimates
    'GPRH',    gprhEst{:};             % 7 columns
    % Row 5: GPRH standard errors
    '',        gprhSE{:};              % 7 columns
    % Row 6: Notes row -- must have 7 cells (pad with extra empty string)
    'Notes:',  'The standard errors are in parentheses.', '*** represents significance level of 1%.', '', '', '', ''
};

% =============== 6) Write to Excel ===============
T = cell2table(table5Data);
writetable(T, 'Table5_GPR_GPRH_Results_v3.xlsx', 'WriteVariableNames', false);

disp('=== Table 5 created: "Table5_GPR_GPRH_Results_v3.xlsx" ===');



%% Table 6: Estimation results of GARCH-MIDAS-GPR model before (2007–2022) 


clc;
clear;

% Pre-conflict subsample (2007–2022)
RTD_pre = xlsread('RTD_GPR_V.xlsx','RTD','C2:C3674');
nobs_pre = size(RTD_pre,1);
EstSample_pre = nobs_pre;
xMonth_pre = xlsread('RTD_GPR_V.xlsx','GPR','C3:C178')./100;
[~, yDate_pre] = xlsread('RTD_GPR_V.xlsx','RTD','A2:A3674');
[~, yDateMonth_pre] = datevec(yDate_pre);
xDay_pre = NaN(nobs_pre,1);
count = 1;
for t = 1:nobs_pre
    if t > 1 && yDateMonth_pre(t) ~= yDateMonth_pre(t-1)
        count = count + 1;
        if count > length(xMonth_pre), break; end
    end
    xDay_pre(t) = xMonth_pre(count);
end

period  = 22;
numLags = 36;
[estParams_pre, EstParamCov_pre] = GarchMidas(RTD_pre, 'Period', period, 'NumLags', numLags, 'X', xDay_pre, 'ThetaM', 1, 'RollWindow', 1, 'estSample', EstSample_pre);
preConflict_params = estParams_pre(:);
preConflict_stdErr = sqrt(diag(EstParamCov_pre));
preCount = length(preConflict_params);

% Post-conflict subsample (2022–2023)
RTD_post = xlsread('RTD_GPR_V.xlsx','RTD','C3675:C4230');
nobs_post = size(RTD_post,1);
EstSample_post = nobs_post;
xMonth_post = xlsread('RTD_GPR_V.xlsx','GPR','C179:C205')./100;
[~, yDate_post] = xlsread('RTD_GPR_V.xlsx','RTD','A3675:A4230');
[~, yDateMonth_post] = datevec(yDate_post);
xDay_post = NaN(nobs_post,1);
count = 1;
for t = 1:nobs_post
    if t > 1 && yDateMonth_post(t) ~= yDateMonth_post(t-1)
        count = count + 1;
        if count > length(xMonth_post), break; end
    end
    xDay_post(t) = xMonth_post(count);
end

[estParams_post, EstParamCov_post] = GarchMidas(RTD_post, 'Period', period, 'NumLags', numLags, 'X', xDay_post, 'ThetaM', 1, 'RollWindow', 1, 'estSample', EstSample_post);
postConflict_params = estParams_post(:);
postConflict_stdErr = sqrt(diag(EstParamCov_post));
postCount = length(postConflict_params);

if preCount ~= postCount
    error('Mismatch: pre = %d, post = %d parameters', preCount, postCount);
end

% Build Table 6
% Define parameter names (modify if necessary)
paramNames = {'\mu','\alpha','\beta','\theta','\omega','m'};
if length(paramNames) ~= preCount
    error('paramNames count (%d) does not match parameter count (%d)', length(paramNames), preCount);
end

formatParams = @(p) sprintf('%.4f', p);
formatStdErr = @(s) sprintf('(%.4f)', s);

numParams = preCount;
preEst = cell(1, numParams);
preSE  = cell(1, numParams);
for i = 1:numParams
    preEst{i} = formatParams(preConflict_params(i));
    preSE{i}  = formatStdErr(preConflict_stdErr(i));
end

postEst = cell(1, numParams);
postSE  = cell(1, numParams);
for i = 1:numParams
    postEst{i} = formatParams(postConflict_params(i));
    postSE{i}  = formatStdErr(postConflict_stdErr(i));
end

% Each row must have 1 + numParams columns = 7 columns.
table6Data = {
    ' ',         paramNames{:};               % Row 1: Header row (7 cells)
    '2007–2022', preEst{:};                   % Row 2: Pre-conflict estimates
    '',          preSE{:};                    % Row 3: Pre-conflict std errors
    '2022–2023', postEst{:};                  % Row 4: Post-conflict estimates
    '',          postSE{:};                   % Row 5: Post-conflict std errors
    'Notes:',    'SE in parentheses', '*, **, *** denote significance at 10%, 5%, 1%', '', '', '', ''   % Row 6: Notes (7 cells)
};

T = cell2table(table6Data);
writetable(T, 'Table6_PrePostConflict_GPR_v2.xlsx', 'WriteVariableNames', false);


%% Table 7: Estimation results of GARCH-MIDAS models with different lag orders

clc;
clear;

% 1) Load data
RTD = xlsread('RTD_GPR_V.xlsx','RTD','C2:C4230');   % daily returns
nobs = size(RTD,1);

% Build daily regressor (xDay) for the GPR model
xMonth = xlsread('RTD_GPR_V.xlsx','GPR','C3:C205')./100;
[~, yDate] = xlsread('RTD_GPR_V.xlsx','RTD','A2:A4230');
[~, yDateMonth] = datevec(yDate);
xDay = NaN(nobs,1);
count = 1;
for t = 1:nobs
    if t > 1 && yDateMonth(t) ~= yDateMonth(t-1)
        count = count + 1;
        if count > length(xMonth)
            break;
        end
    end
    xDay(t) = xMonth(count);
end

% 2) Define lag orders and enforce two omega terms
kValues = [12, 24, 36, 48];
beta2Flag = 1;  % Enforce two omega terms

% 3) Set parameter names (7 total)
% Desired order: mu, alpha, beta, theta, omega_1, omega_2, m
paramNames = {'\mu','\alpha','\beta','\theta','\omega_1','\omega_2','m'};
numParams  = length(paramNames);  % should be 7

% 4) Formatting functions
formatParams = @(p) sprintf('%.4f', p);
formatStdErr = @(s) sprintf('(%.4f)', s);

% 5) Initialize cell array for the final table
table7Data = {};

% 6) Loop over each lag order k
for k = kValues
    % --- 6.1) Estimate RV model ---
    [estParams_RV, EstParamCov_RV, ~, ~] = GarchMidas(...
        RTD, 'Period', 22, 'NumLags', k, 'ThetaM', 1, 'RollWindow', 1, ...
        'estSample', 4229, 'Beta2Para', beta2Flag);
    rvParams = estParams_RV(:);
    rvStdErr = sqrt(diag(EstParamCov_RV));
    if length(rvParams) ~= numParams
        error('RV model with k=%d returned %d params, expected %d.', k, length(rvParams), numParams);
    end
    
    % --- 6.2) Estimate GPR model ---
    [estParams_GPR, EstParamCov_GPR, ~, ~] = GarchMidas(...
        RTD, 'Period', 22, 'NumLags', k, 'X', xDay, 'ThetaM', 1, 'RollWindow', 1, ...
        'estSample', 4229, 'Beta2Para', beta2Flag);
    gprParams = estParams_GPR(:);
    gprStdErr = sqrt(diag(EstParamCov_GPR));
    if length(gprParams) ~= numParams
        error('GPR model with k=%d returned %d params, expected %d.', k, length(gprParams), numParams);
    end
    
    % --- 6.3) Estimate RV+GPR model ---
    x21 = [RTD, xDay];
    [estParams_RVGPR, EstParamCov_RVGPR, ~, ~] = GarchMidas2(...
        RTD, 'X', x21, 'Period', 22, 'NumLags', k, 'LogTau', 1, 'EstSample', 4229, ...
        'RollWindow', 1, 'thetaM', 1, 'Beta2Para', beta2Flag);
    % If extra parameters are returned, keep only the first 7.
    if length(estParams_RVGPR) > numParams
        estParams_RVGPR = estParams_RVGPR(1:numParams);
        EstParamCov_RVGPR = EstParamCov_RVGPR(1:numParams, 1:numParams);
    end
    rvgprParams = estParams_RVGPR(:);
    rvgprStdErr = sqrt(diag(EstParamCov_RVGPR));
    if length(rvgprParams) ~= numParams
        error('RV+GPR model with k=%d returned %d params, expected %d.', k, length(rvgprParams), numParams);
    end
    
    % --- 6.4) Convert numerical results to cell strings ---
    rvEst = cell(1, numParams); rvSE = cell(1, numParams);
    gprEst = cell(1, numParams); gprSE = cell(1, numParams);
    rgEst = cell(1, numParams); rgSE = cell(1, numParams);
    for i = 1:numParams
        rvEst{i} = formatParams(rvParams(i));
        rvSE{i}  = formatStdErr(rvStdErr(i));
        gprEst{i} = formatParams(gprParams(i));
        gprSE{i}  = formatStdErr(gprStdErr(i));
        rgEst{i} = formatParams(rvgprParams(i));
        rgSE{i}  = formatStdErr(rvgprStdErr(i));
    end
    
    % --- 6.5) Build block for current lag order (each block: 7 rows, 8 columns) ---
    % Use curly braces to create a cell array row.
    block = {
        ['k=', num2str(k)], paramNames{:};         % Row 1: Header row
        'RV',     rvEst{:};                        % Row 2: RV estimates
        '',       rvSE{:};                         % Row 3: RV std errors
        'GPR',    gprEst{:};                       % Row 4: GPR estimates
        '',       gprSE{:};                        % Row 5: GPR std errors
        'RV+GPR', rgEst{:};                        % Row 6: RV+GPR estimates
        '',       rgSE{:}                          % Row 7: RV+GPR std errors
    };
    
    % Append this block to table7Data
    table7Data = [table7Data; block];  %#ok<AGROW>
end

% 7) Add a final notes row (8 cells)
notesRow = {'Notes:', 'k=12,24,36,48 are the lag orders.', 'Two omega terms enforced.', '***, **, * indicate significance at 1%, 5%, 10%.', '', '', '', ''};
table7Data = [table7Data; notesRow];

% 8) Write to Excel (without auto-generated headers)
T = cell2table(table7Data);
writetable(T, 'Table7_LagOrders_Results.xlsx', 'WriteVariableNames', false);

disp('=== Table 7 created: "Table7_LagOrders_Results.xlsx" with 2 omega terms ===');

clc;
clear;

%% Table 8: Estimation results of GARCH-MIDAS models for oil-exporting and 

clc;
clear;

% 1) Common Settings
RTD = xlsread('RTD_GPR_V.xlsx','RTD','C2:C4230');  % daily returns
nobs = size(RTD,1);

period    = 22;
numLags   = 36;
Beta2Para = 0;   % => 1 "omega" parameter, total 6 parameters
thetaM    = 1;
Rollwindow= 1;
estSample = 4229;

paramNames = {'\mu','\alpha','\beta','\theta','\omega','m'};
numParams  = length(paramNames);  % 6

formatParams = @(p) sprintf('%.4f', p);
formatStdErr = @(s) sprintf('(%.4f)', s);

% 2) Build daily series for each country
% We have 6 columns in "GPR_Country": B->G
%  - B3:B205 = xMonth_CAN
%  - C3:C205 = xMonth_RUS
%  - D3:D205 = xMonth_BRA
%  - E3:E205 = xMonth_USA
%  - F3:F205 = xMonth_JPN
%  - G3:G205 = xMonth_IDN

xMonth_CAN = xlsread('RTD_GPR_V.xlsx','GPR_Country','B3:B205')./100;
xMonth_RUS = xlsread('RTD_GPR_V.xlsx','GPR_Country','C3:C205')./100;
xMonth_BRA = xlsread('RTD_GPR_V.xlsx','GPR_Country','D3:D205')./100;
xMonth_USA = xlsread('RTD_GPR_V.xlsx','GPR_Country','E3:E205')./100;
xMonth_JPN = xlsread('RTD_GPR_V.xlsx','GPR_Country','F3:F205')./100;
xMonth_IDN = xlsread('RTD_GPR_V.xlsx','GPR_Country','G3:G205')./100;

% We'll build a helper function to expand monthly xMonth -> daily xDay
[~, yDate] = xlsread('RTD_GPR_V.xlsx','RTD','A2:A4230');
[~, yDateMonth] = datevec(yDate);

buildDaily = @(xMonthData) expandMonthlyToDaily(xMonthData, yDateMonth);

% Now build xDay for each country
xDay_CAN = buildDaily(xMonth_CAN);
xDay_RUS = buildDaily(xMonth_RUS);
xDay_BRA = buildDaily(xMonth_BRA);
xDay_USA = buildDaily(xMonth_USA);
xDay_JPN = buildDaily(xMonth_JPN);
xDay_IDN = buildDaily(xMonth_IDN);

% 3) Estimate each country's GARCH-MIDAS
% 3.1) Oil-exporting: CAN, RUS, BRA
[canParams, canSE] = estimateGarchMidasCountry(RTD, xDay_CAN, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM);
[rusParams, rusSE] = estimateGarchMidasCountry(RTD, xDay_RUS, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM);
[braParams, braSE] = estimateGarchMidasCountry(RTD, xDay_BRA, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM);

% 3.2) Oil-importing: USA, JPN, IDN
[usaParams, usaSE] = estimateGarchMidasCountry(RTD, xDay_USA, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM);
[jpnParams, jpnSE] = estimateGarchMidasCountry(RTD, xDay_JPN, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM);
[idnParams, idnSE] = estimateGarchMidasCountry(RTD, xDay_IDN, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM);

% 4) Convert numeric -> cell arrays of strings
[canEst, canStd] = convertParamsSE(canParams, canSE, formatParams, formatStdErr);
[rusEst, rusStd] = convertParamsSE(rusParams, rusSE, formatParams, formatStdErr);
[braEst, braStd] = convertParamsSE(braParams, braSE, formatParams, formatStdErr);

[usaEst, usaStd] = convertParamsSE(usaParams, usaSE, formatParams, formatStdErr);
[jpnEst, jpnStd] = convertParamsSE(jpnParams, jpnSE, formatParams, formatStdErr);
[idnEst, idnStd] = convertParamsSE(idnParams, idnSE, formatParams, formatStdErr);

% 5) Build Table 8 as a single cell array
% Each row has 1 (label) + 6 (params) = 7 cells
table8Data = {
    % Panel A heading
    'Panel A: Oil-exporting countries', paramNames{:};
    % CAN
    'CAN', canEst{:};
    '',    canStd{:};
    % RUS
    'RUS', rusEst{:};
    '',    rusStd{:};
    % BRA
    'BRA', braEst{:};
    '',    braStd{:};
    
    % Panel B heading
    'Panel B: Oil-importing countries', paramNames{:};
    % USA
    'USA', usaEst{:};
    '',    usaStd{:};
    % JPN
    'JPN', jpnEst{:};
    '',    jpnStd{:};
    % IDN
    'IDN', idnEst{:};
    '',    idnStd{:};
};

% 6) Write to Excel
T = cell2table(table8Data);
writetable(T, 'Table8_OilCountries_Results.xlsx', 'WriteVariableNames', false);

disp('=== Table 8 created: "Table8_OilCountries_Results.xlsx" ===');

%            Helper Functions

function xDay = expandMonthlyToDaily(xMonthData, yDateMonth)
    % Expand monthly series xMonthData -> daily xDay
    nobs = length(yDateMonth);
    xDay = NaN(nobs,1);
    count = 1;
    xLen = length(xMonthData);
    for ii = 1:nobs
        if ii>1 && yDateMonth(ii) ~= yDateMonth(ii-1)
            count = count + 1;
            if count> xLen, break; end
        end
        xDay(ii) = xMonthData(count);
    end
end

function [paramsVec, stdErrVec] = estimateGarchMidasCountry(RTD, xDay, ...
    period, numLags, Beta2Para, estSample, Rollwindow, thetaM)
    % Single-country GarchMidas call returning 6 parameters + std. errors
    [estParams, EstParamCov, ~, ~] = GarchMidas( ...
        RTD, ...
        'Period',       period, ...
        'NumLags',      numLags, ...
        'X',            xDay, ...
        'ThetaM',       thetaM, ...
        'RollWindow',   Rollwindow, ...
        'estSample',    estSample, ...
        'Beta2Para',    Beta2Para ...
    );
    paramsVec = estParams(:);
    stdErrVec = sqrt(diag(EstParamCov));
    if length(paramsVec)~=6
        error('Model returned %d parameters, expected 6. Check Beta2Para.', length(paramsVec));
    end
end

function [estCell, stdCell] = convertParamsSE(paramsVec, stdErrVec, fEst, fStd)
    % Convert numeric param & std error vectors -> cell arrays of strings
    n = length(paramsVec);
    estCell = cell(1,n);
    stdCell = cell(1,n);
    for i = 1:n
        estCell{i} = fEst(paramsVec(i));
        stdCell{i} = fStd(stdErrVec(i));
    end
end

%% Table 9: Estimation results of GARCH-MIDAS models in two crisis periods:

% =============== 1) First crisis period: Aug 2013–Apr 2018 ===============
% Subsample of RTD
RTD1 = xlsread('RTD_GPR_V.xlsx','RTD','C1566:C2748');
nobs1 = size(RTD1,1);
EstSample1 = nobs1;  % in-sample only

% Build xDay1 from GPR sheet
xMonth1 = xlsread('RTD_GPR_V.xlsx','GPR','C78:C134')./100;
[~, yDate1] = xlsread('RTD_GPR_V.xlsx','RTD','A1566:A2748');
[~, yDateMonth1] = datevec(yDate1);
xDay1 = expandMonthlyToDaily(xMonth1, yDateMonth1);

% Estimate
[estParams1, EstParamCov1, ~, ~] = GarchMidas( ...
    RTD1, ...
    'Period',     22, ...
    'NumLags',    4,  ...  % as in your code
    'X',          xDay1, ...
    'ThetaM',     1, ...
    'RollWindow', 1, ...
    'estSample',  EstSample1, ...
    'Beta2Para',  0 ...
);
params1 = estParams1(:);
stdErr1 = sqrt(diag(EstParamCov1));

% =============== 2) Second crisis period: Jan 2020–Oct 2021 ===============
RTD2 = xlsread('RTD_GPR_V.xlsx','RTD','C3171:C3631');
nobs2 = size(RTD2,1);
EstSample2 = nobs2;

% Build xDay2 from GPR sheet
xMonth2 = xlsread('RTD_GPR_V.xlsx','GPR','C155:C176')./100;
[~, yDate2] = xlsread('RTD_GPR_V.xlsx','RTD','A3171:A3631');
[~, yDateMonth2] = datevec(yDate2);
xDay2 = expandMonthlyToDaily(xMonth2, yDateMonth2);

% Estimate
[estParams2, EstParamCov2, ~, ~] = GarchMidas( ...
    RTD2, ...
    'Period',     22, ...
    'NumLags',    3, ...  % as in your code
    'X',          xDay2, ...
    'ThetaM',     1, ...
    'RollWindow', 1, ...
    'estSample',  EstSample2, ...
    'Beta2Para',  0 ...
);
params2 = estParams2(:);
stdErr2 = sqrt(diag(EstParamCov2));

% =============== 3) Convert numeric -> strings and build Table 9 ===============
% We assume the model has 6 parameters:
paramNames = {'\mu','\alpha','\beta','\theta','\omega','m'};
numParams  = length(paramNames);

formatEst = @(p) sprintf('%.4f', p);
formatSE  = @(s) sprintf('(%.4f)', s);

% 3.1) Check we have 6 parameters in each model
if length(params1) ~= numParams
    error('First crisis period returned %d params, expected %d.', length(params1), numParams);
end
if length(params2) ~= numParams
    error('Second crisis period returned %d params, expected %d.', length(params2), numParams);
end

% 3.2) Build row by row
% Row 1: blank in col 1, then the two column headers
table9Data = {
    ' ', 'August 2013–April 2018', 'January 2020–October 2021'
};

% Then for each parameter i, we add 2 rows:
%   Row (2i): paramName, param1, param2
%   Row (2i+1): blank, stdErr1, stdErr2
for i = 1:numParams
    thisParamName = paramNames{i};
    est1Str = formatEst(params1(i));
    est2Str = formatEst(params2(i));
    se1Str  = formatSE(stdErr1(i));
    se2Str  = formatSE(stdErr2(i));
    
    % Row for estimates
    table9Data = [ table9Data;
        { thisParamName, est1Str, est2Str } ];
    % Row for std. errors
    table9Data = [ table9Data;
        { '', se1Str, se2Str } ];
end

% Finally, add a notes row
notesRow = {
    'Notes:', 'The standard errors are in parentheses. * and *** denote significance levels of 10% and 1%.', ''
};
table9Data = [table9Data; notesRow];

% =============== 4) Write to Excel ===============
T = cell2table(table9Data);
writetable(T, 'Table9_CrisisPeriods_Results.xlsx', 'WriteVariableNames', false);

disp('=== Table 9 created: "Table9_CrisisPeriods_Results.xlsx" ===');



% Helper function to expand monthly -> daily

function xDay = expandMonthlyToDaily(xMonthData, yDateMonth)
    nobs = length(yDateMonth);
    xDay = NaN(nobs,1);
    count = 1;
    for ii = 1:nobs
        if ii>1 && yDateMonth(ii)~=yDateMonth(ii-1)
            count = count + 1;
            if count>length(xMonthData)
                break;
            end
        end
        xDay(ii) = xMonthData(count);
    end
end
