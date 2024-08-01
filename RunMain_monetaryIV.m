clear; clc;
addpath([pwd, '\Functions'])
%% Settings
nlags = 6; % VAR lags
irfperiods = 48; % periods to generate IRFs for
nsave = 1000; % Final number of Bayes draws to save
nburn = 100; % Burn in   

%% Import data to estimate a VAR containing: log GDP, Log core CPI, unemployment, fed funds, log c&icredit, sloos C&I standards, AAAcorp spread

% Estimate with shadow rate
[USDataQ, meta]=xlsread('MonthlydataMonPol.xlsx','Import'); 

% Cutoff point for estimation of the VAR (decomp will be dont on whole dataset)
loceststart = find(strcmp(meta(2:end,1),'1/31/1978'));
locest = find(strcmp(meta(2:end,1),'12/31/2007'));

IVGK = USDataQ(loceststart+nlags:locest,end-1);
IVMAR = USDataQ(loceststart+nlags:locest,end);
USDataQ = USDataQ(loceststart:locest,1:end-2);

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

%% 
idents = {'IVGK', 'IVMAR'};

%% Estimate the BVAR
Ident = 7; % Ident 8 in "BootVAR" is the IV identification. It requires that the VAR variable ordered first is the one being instrumented
options = [];


% Estimate Gertler & Karadi IV
options.IV = IVGK;
[InvA_draws.IVGK ,ALPHA_draws.IVGK, SIGMA_draws.IVGK, HDshock_draws.IVGK, HDinit_draws.IVGK, HDconst_draws.IVGK, IRF_draws.IVGK,FEVD_draws.IVGK] = ...
    BootVAR(X, Y, nlags, nvars, Ident, nsave, irfperiods, options);

% Estimate Miranda-Aggrippino & Ricco IV
options.IV = IVMAR;
[InvA_draws.IVMAR ,ALPHA_draws.IVMAR, SIGMA_draws.IVMAR, HDshock_draws.IVMAR, HDinit_draws.IVMAR, HDconst_draws.IVMAR, IRF_draws.IVMAR,FEVD_draws.IVMAR] = ...
    BootVAR(X, Y, nlags, nvars, Ident, nsave, irfperiods, options);


%% Plot IRFs
FillColor   =[.85 .85 .85];
figure
FontSizeset = 10;

for ii = 1:size(idents,2) 
subplot(size(Y,2),size(idents,2),ii*1)
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(idents{ii})(:,1,:,1)),[0.16])  fliplr(quantile(squeeze(IRF_draws.(idents{ii})(:,1,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.(idents{ii})(:,1,:,1)),1)','-k','LineWidth',0.5)
    hold on
   plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    text = [idents{ii}];
    title(text)
    ylabel(['Fed funds'])
    set(gca,'FontSize', FontSizeset)
 subplot(size(Y,2),size(idents,2),(size(idents,2)+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ cumsum(quantile(squeeze(IRF_draws.(idents{ii})(:,2,:,1)),[0.16]))  fliplr(cumsum(quantile(squeeze(IRF_draws.(idents{ii})(:,2,:,1)),[0.84])))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(cumsum(nanmedian(squeeze(IRF_draws.(idents{ii})(:,2,:,1)),1)'),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Core CPI'])
    set(gca,'FontSize', FontSizeset)
  subplot(size(Y,2),size(idents,2),(size(idents,2)*2+ii))
   fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(idents{ii})(:,3,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(idents{ii})(:,3,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.(idents{ii})(:,3,:,1)),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['IP'])
    set(gca,'FontSize', FontSizeset)
   subplot(size(Y,2),size(idents,2),(size(idents,2)*3+ii))
      fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(idents{ii})(:,4,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(idents{ii})(:,4,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.(idents{ii})(:,4,:,1)),1),'-k','LineWidth',0.5)
    hold on

    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Consumption'])
    
    set(gca,'FontSize', FontSizeset)
       subplot(size(Y,2),size(idents,2),(size(idents,2)*4+ii))
      fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(idents{ii})(:,5,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(idents{ii})(:,5,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.(idents{ii})(:,5,:,1)),1),'-k','LineWidth',0.5)
    hold on

    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Unemp'])
    
    set(gca,'FontSize', FontSizeset)

           subplot(size(Y,2),size(idents,2),(size(idents,2)*5+ii))
      fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(idents{ii})(:,6,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(idents{ii})(:,6,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.(idents{ii})(:,6,:,1)),1),'-k','LineWidth',0.5)
    hold on

    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['C&I Credit'])
    
    set(gca,'FontSize', FontSizeset)
               subplot(size(Y,2),size(idents,2),(size(idents,2)*6+ii))
      fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(idents{ii})(:,7,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(idents{ii})(:,7,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(nanmedian(squeeze(IRF_draws.(idents{ii})(:,7,:,1)),1),'-k','LineWidth',0.5)
    hold on

    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Corp Spreads'])
    
    set(gca,'FontSize', FontSizeset)
end

set(gcf, 'Position',  [100, 100, 800, 800])

