clear; clc;
addpath([pwd, '/Functions'])
%% Settings

irfperiods = 48; % periods to generate IRFs for
nsave = 1000; % Final number of Bayes draws to save
nburn = 100; % Burn in   

%% Import data to estimate a VAR containing: log GDP, Log core CPI, unemployment, fed funds, log c&icredit, sloos C&I standards, AAAcorp spread

[USDataQ]=readtable('MonthlydataMonPol.csv');


loceststart = find(USDataQ.Var1==197901);
locest = find(USDataQ.Var1==201412);

USDataQ = USDataQ{:,2:end};

IVGK = USDataQ(loceststart:locest,end-1);
IVMAR = USDataQ(loceststart:locest,end);
USDataQ = USDataQ(loceststart:locest,1:end-2);

%% Use AIC to choose lag length
maxlag = 12;
nlags = AICchooselag(USDataQ,maxlag);

%% Create data
[T,nvars]=size(USDataQ);
X = [];

for p=1:nlags
    X(:,1+(p-1)*nvars:p*nvars)=USDataQ((nlags+1-p):(T-p),:);
    %gives (T-nlags) x nvars*nlags matrix of lags for all variables                                         
end
% Add constant
const = ones(T-nlags,1);
nobs = T-nlags;
X = [X const];
%rescaling Y variables since we loose nlags observations
Y=USDataQ((nlags+1):end,:);
IVMAR = IVMAR(nlags+1:end,1);
IVGK = IVGK(nlags+1:end,1);
%% 
idents = {'IVGK', 'IVMAR'};

%% Weak IV tests - note that this is instrumenting the policy rate rather than FF4 as in the original papers
disp("Miranda-Agrippino-Ricco:")
Fstat = WeakIVTest(Y(1:end,1),X, IVMAR);
disp("Gertler-Karadi:")
Fstat = WeakIVTest(Y(1:end,1),X, IVGK);


%% Estimate the BVAR
Ident = 7; % Ident 8 in "BootVAR" is the IV identification. It requires that the VAR variable ordered first is the one being instrumented
options = [];


% Estimate Gertler & Karadi IV
options.IV = IVGK;
[InvA_draws.IVGK ,ALPHA_draws.IVGK, SIGMA_draws.IVGK, HDshock_draws.IVGK, HDinit_draws.IVGK, HDconst_draws.IVGK, IRF_draws.IVGK,FEVD_draws.IVGK] = ...
    BootVAR(X, Y, nlags, nvars, Ident, nsave, irfperiods, options);

% Estimate Miranda-Agrippino & Ricco IV
options.IV = IVMAR;
[InvA_draws.IVMAR ,ALPHA_draws.IVMAR, SIGMA_draws.IVMAR, HDshock_draws.IVMAR, HDinit_draws.IVMAR, HDconst_draws.IVMAR, IRF_draws.IVMAR,FEVD_draws.IVMAR] = ...
    BootVAR(X, Y, nlags, nvars, Ident, nsave, irfperiods, options);



%% Plot IRFs
% Names to use in use for subplot chart titles
varnames = {'Fed Funds', 'Core CPI', 'IP', 'Consumption', 'Unemployment', 'Credit', 'Corp spreads'};
% Order and selection of vars to plot IRFs
varuse = [1,2,3,4,5,6,7];
% Vars measured in log-differences that will be cumulated for the plot to
% measure level impact (Core CPI in this case)
diffs = [0,1,0,0,0,0,0];

% Shock to plot 
shock = 1; % IV shock is ordered first

%Title for whole figure
wholetitle = "Miranda-Agrippino-Ricco IV";
PlotSVARSingle(IRF_draws.IVMAR, varnames, varuse, diffs, irfperiods, wholetitle, shock)

set(gcf, 'Position',  [100, 100, 800, 800])

wholetitle = "Gertler-Karadi IV";
PlotSVARSingle(IRF_draws.IVGK, varnames, varuse, diffs, irfperiods, wholetitle, shock)
set(gcf, 'Position',  [100, 100, 800, 800])
