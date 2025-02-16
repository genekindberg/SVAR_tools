function Fstat = WeakIVTest(Y, X, Z)

% Y is endogenous variable of interest, X lags of all exogenous variables, Z is instrument
% remove nan rows (often present given IV will not cover whole period) 
t = any(isnan([Y, X,Z]),2);
Y(t,:) = [];
X(t,:) = [];
Z(t,:) = [];


% First stage
Zstage1 = [Z X]; % Including instruments and endogenous


Beta = Zstage1\Y; % Coefficients from the first stage
Y_hat = Zstage1 * Beta; % Predicted values of X (X_hat)

% Calculate the residuals from the first-stage regression
residuals_first_stage = Y - Y_hat;
% Calculate the SSR for the unrestricted model (with instruments)
SSR_unrestricted = sum(residuals_first_stage.^2);


% Calculate the Sum of Squared Residuals (SSR) for the restricted model (without instruments).
Zstage2 = [X]; % Including instruments and constant
Beta = Zstage2\Y; % Coefficients from the first stage
Y_hatstage2 = Zstage2 * Beta; % Predicted values of X (X_hat)
SSR_restricted = sum((Y - Y_hatstage2).^2);




% Number of instruments (excluding the constant)
k = 1;

% Number of observations
n = size(X,1);

% Number of endogenous variables (including the constant)
m = size(X,2);

% First-stage F-statistic
Fstat = ((SSR_restricted - SSR_unrestricted) / k) / (SSR_unrestricted / (n - m - k));

X = sprintf('F statistic: %.1f',Fstat);
disp(X)
end