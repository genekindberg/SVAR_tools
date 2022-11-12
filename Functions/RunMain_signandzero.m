clear; clc;
addpath([pwd, '\Functions'])
% Script estimates a technology shock, demand shock, monetary
% policy shock, and supply shock using sign and zero restrictions, and
% a restriction that demand shocks have a larger share of the FEVD of labor
% productivity in the first year than technology, and a smaller share at 5
% years.
%% Settings
nlags = 4; % VAR lags
irfperiods = 20; % periods to generate IRFs for
nsave = 1000; % Final number of Bayes draws to save
nburn = 500; % Burn in   

%% Import data to estimate a VAR containing: log labor productivity, log hours worked, log investment/GDP, log consumption/GDP, consumer price inflation, long-term treasury yields.
USDataQ=xlsread('SVARs.xls','Quarterly');
USDataM=xlsread('SVARs.xls','Monthly');
USDataQ = USDataQ(5:end,:); % get rid of 1947 to be consistent with monthly data

% Manipulations for data
% Monthly to quarterly
CivEmp = mean([USDataM(1:3:end-2,1) USDataM(2:3:end-1,1) USDataM(3:3:end,1)],2); %Civilian employment -CE16OV
CivPop = mean([USDataM(1:3:end-2,2) USDataM(2:3:end-1,2) USDataM(3:3:end,2)],2); % Civilian population - CNP16OV
LRates = mean([USDataM(1:3:end-2,3) USDataM(2:3:end-1,3) USDataM(3:3:end,3)],2); % Civilian population - 10-year treasuries - GS10 1953 Q2 onwards

% Quarterly
RGDPpc = USDataQ(:,1); % GDP per capita - A939RX0Q048SBEA
AWH = USDataQ(:,9); % Average weekly hours - PRS85006023
GDP = USDataQ(:,3); % Nominal GDP - GDP
PCEDef = USDataQ(:,2); % PCE deflator - DPCERD3Q086SBEA 
Lprod = USDataQ(:,5); % Labor productivity


% Manipulations
H=log(AWH.*CivEmp./CivPop)*100; % log hours worked
Lprod = log(Lprod)*100; % log labor productivity levels
PI = (log(PCEDef(2:end,1))-log(PCEDef(1:end-1,1)))*100; % log inflation

% log-differenced labor productivity and hours worked, quarterly inflation.
% and 10-year bond yields.
DatasetQdiffs = [(Lprod(2:end,1)-Lprod(1:end-1,1)) (H(2:end,1)-H(1:end-1,1)) PI(1:end,1) LRates(2:end,1)];

% Start from 1953 Q2:2018Q3
DatasetQdiffs = DatasetQdiffs(21:end-1,:); 



%% Create data
[T,nvars]=size(DatasetQdiffs);
X = [];

for p=1:nlags
    X(:,1+(p-1)*nvars:p*nvars)=DatasetQdiffs((nlags+1-p):(T-p),:);
    %gives (T-nlags) x nvars*nlags matrix of lags for all variables                                         
end
const = ones(T-nlags,1);
nobs = T-nlags;
X = [X const];
%rescaling variables since we loose nlags observations
Y=DatasetQdiffs((nlags+1):end,:);

%% Set sign restrictions
shocks = {'Technology', 'Demand', 'Mon. Pol.', 'Supply'};
% Sign restrictions for each shock (Shocks 1-4)
% [Shock, Variable, Period, Pos(1)/Neg(-1), Magnitude (at least x response)]

restspec = [1 1 inf 1 0; % Tech shock has positive long run impact
                          2 2 0 1 0;% Demand boosts employment
                          2 3 0 1 0.2;% Demand boosts inflation
                          2 4 0 1 0; % Demand boosts interest
                          3 2 0 -1 0;% MonPol reduces employment
                          3 3 0 -1 0;% MonPol reduces inflation
                          3 4 0 1 0; % MonPol boosts interest
                          4 2 0 -1 0;% Supply reduces employment
                          4 3 0 1 0;% Supply boost inflation
                          4 4 0 1 0]; % Supply boosts interest] 

Exclusion = [2 1 inf; % Demand no long run impact
             3 1 inf; % Mon Pol no long run impact
             4 1 inf]; % Supply no long run impact
% Exclusion = [];

%% Impose ratio restriction
% [Shock, variable 1, variable 2, horizon, sign, magnitude] so that variable
% 1/variable 2 >(1), <(-1) magnitude
% For example a monetary policy shock that has the opposite impact on
% inflation and interest rates
% ratio = [3 4 3 0 -1 0; % ratio of impact on rates over impact on
% inflation is less than zero.
%           ]; 
ratio = [];

%% Assign to function arguments
S = restspec;
E = Exclusion;


%% Set FEVD restrictions
% % demand shock has larger impact on labor productivity FEVD in first year, technology at 5 years
% FEVDperiods = [4,20]; % 1 year and 5 years
% FEVDrestrict = {[1 2 1],[1 1 2]}; % {period 1(endogenous variable, shock that is biggest, shock that is smallest FEVD), period 2(endogenous variable, shock that is biggest, shock that is smallest FEVD)}
% % multiple restrictions in each period can be stacked i.e. [1 2 1; 1 3 1]
% % would impose that the second and third shocks had a larger share of FEVD
% % for variable 1 than the first shock in period 1.
FEVDperiods = [];
FEVDrestrict = [];
%% Estimate the BVAR
Ident = 4; % the fourth identification is the one that allows sign, zero and FEVD restrictions.
options = [];
options.restrictsign = S; % add set of sign restrictions
options.restrictzero = E; % add set of zero restrictions
options.trys = 1000; % how many random matrices are drawn to try and fit the sign and zero restrictions for each bayesian draw of coefficients and variance before giving up and redrawing coefficients.
options.FEVDperiods = FEVDperiods;
options.FEVDrestrict = FEVDrestrict;
options.ratio = ratio;

% BVAR priors
lambdas = zeros(4,1);
lambdas(1) = 1;% tightness on own AR1 lags
lambdas(2) = 1;% tightness on lags of other variables
lambdas(3) = 1;% tightness on own additional lags - how much tighter prior gets over time
lambdas(4) = 10^2; % tightness on constant and exogenous

[InvA_draws.Long_RunMixed ,ALPHA_draws.Long_RunMixed, SIGMA_draws.Long_RunMixed, HDshock_draws.Long_RunMixed, HDinit_draws.Long_RunMixed, HDconst_draws.Long_RunMixed, IRF_draws.Long_RunMixed,FEVD_draws.Long_RunMixed] = ...
    BVAR(X, Y, nlags, nvars, Ident, nburn, nsave, irfperiods, lambdas, options);


%% Plot
FillColor   =[.85 .85 .85];
figure


for ii = 1:size(shocks,2) 
subplot(size(Y,2),size(shocks,2),ii*1)
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(cumsum(squeeze(IRF_draws.Long_RunMixed(:,1,:,ii))')',[0.16])  fliplr(quantile(cumsum(squeeze(IRF_draws.Long_RunMixed(:,1,:,ii))')',[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(cumsum(squeeze(IRF_draws.Long_RunMixed(:,1,:,ii))'),2),'-k','LineWidth',0.5)
    hold on
   plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    text = [shocks{ii}];
    title(text)
    ylabel(['Labour Productivity'])
 subplot(size(Y,2),size(shocks,2),(size(shocks,2)+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(cumsum(squeeze(IRF_draws.Long_RunMixed(:,2,:,ii))')',[0.16])  fliplr(quantile(cumsum(squeeze(IRF_draws.Long_RunMixed(:,2,:,ii))')',[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(cumsum(squeeze(IRF_draws.Long_RunMixed(:,2,:,ii))'),2),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Hours'])
  subplot(size(Y,2),size(shocks,2),(size(shocks,2)*2+ii))
   fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.Long_RunMixed(:,3,:,ii)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.Long_RunMixed(:,3,:,ii)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.Long_RunMixed(:,3,:,ii)),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Inflation'])
   subplot(size(Y,2),size(shocks,2),(size(shocks,2)*3+ii))
      fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.Long_RunMixed(:,4,:,ii)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.Long_RunMixed(:,4,:,ii)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.Long_RunMixed(:,4,:,ii)),1),'-k','LineWidth',0.5)
    hold on

    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Long Rates'])
end

set(gcf, 'Position',  [100, 100, 1000, 800])

