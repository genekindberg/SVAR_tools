clear; clc;
addpath([pwd, '\Functions'])
% Script estimates a several technology shock identifications - then
% produces a range of outputs, including IRFs, FEVD, historical
% decompositions. Further identifications will be added in due course.
%% Settings
nlags = 4; % VAR lags
irfperiods = 40; % periods to generate IRFs for
nsave = 1000; % Final number of draws to save
nburn = 500; % Burn in   

% VARs being estimated - using for chart loops - can be expanded
VAR_Or = {'Long_Run','Max_Share'};
VAR_Nice = {'Long-Run','Max Share'}; % text for printing on charts
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
PCDG = USDataQ(:,6); % Durable goods consumption PCDG
GPDI = USDataQ(:,4); % Private investment = GPDI
PCND = USDataQ(:,8); % Non-durable PCND
PCESV = USDataQ(:,7); %  Services expenditure (consumption) PCESV
GDP = USDataQ(:,3); % Nominal GDP - GDP
PCEDef = USDataQ(:,2); % PCE deflator - DPCERD3Q086SBEA 
Lprod = USDataQ(:,5); % Labor productivity
% Manipulations
I_Y = log(100*((PCDG+GPDI)./GDP))*100;
C_Y = log(100*((PCND+PCESV)./GDP))*100;

% Manipulations
H=log(AWH.*CivEmp./CivPop)*100;
Lprod = log(Lprod)*100;
PI = (log(PCEDef(2:end,1))-log(PCEDef(1:end-1,1)))*100;

DatasetQ = [Lprod(2:end,1) H(2:end,1) I_Y(2:end,1) C_Y(2:end,1) PI(1:end,1) LRates(2:end,1)];
DatasetQGali = [(Lprod(2:end,1)-Lprod(1:end-1,1)) (H(2:end,1)-H(1:end-1,1)) I_Y(2:end,1) C_Y(2:end,1) PI(1:end,1) LRates(2:end,1)];

% Start from 1953 Q2:2018Q3
DatasetQ = DatasetQ(21:end-1,:); % Productivity data in levels
DatasetQGali = DatasetQGali(21:end-1,:); % Gali uses differenced productivity data


%% Make data for VARs

X = [];
[T,nvars]=size(DatasetQ);
for p=1:nlags
    X(:,1+(p-1)*nvars:p*nvars)=DatasetQ((nlags+1-p):(T-p),:); %gives (T-nlags) x nvars*nlags matrix of lags for all variables                                         
end

const = ones(T-nlags,1);
X = [X const]; % Constant always goes at the end
Y=DatasetQ((nlags+1):T,:);% rescaling variables since we loose nlags observations through the lagging 
nobs = size(Y,1);

XGali = [];
[T,nvars]=size(DatasetQ);
for p=1:nlags
    XGali(:,1+(p-1)*nvars:p*nvars)=DatasetQGali((nlags+1-p):(T-p),:);
    %gives (T-nlags) x nvars*nlags matrix of lags for all variables                                         
end

const = ones(T-nlags,1);
XGali = [XGali const];
%rescaling variables since we loose nlags observations through the lagging 
YGali=DatasetQGali((nlags+1):T,:);
nobs = size(YGali,1);


%% Estimate VARs
optionsLR = [];
IdentLR = 1; % Long-run restriction - always establishes shock that leads to long-run impact on first variable
optionsLR.diff=1; % Cumulates FEVD so that it is expressed as contribution to variance of variable in levels
optionsLR.target=1; % So function knows which FEVD to cumulate - productivity in differences.
[InvA_draws.Long_Run ,ALPHA_draws.Long_Run, SIGMA_draws.Long_Run, HDshock_draws.Long_Run, HDinit_draws.Long_Run, HDconst_draws.Long_Run, IRF_draws.Long_Run,FEVD_draws.Long_Run] = ...
    BVAR(XGali, YGali, nlags, nvars, IdentLR, nburn, nsave, irfperiods, optionsLR);

IdentMS = 2; % Max-Share restriction
optionsMS = [];
optionsMS.target = 1; % which var is target of maximization
optionsMS.horz = 40; % horizon for max FEVD 40 quarters
optionsMS.diff=0;
optionsMS.rotatesign=1; % Flip identified shocks to be positive
[InvA_draws.Max_Share ,ALPHA_draws.Max_Share, SIGMA_draws.Max_Share, HDshock_draws.Max_Share, HDinit_draws.Max_Share, HDconst_draws.Max_Share, IRF_draws.Max_Share,FEVD_draws.Max_Share] = ...
    BVAR(X, Y, nlags, nvars, IdentMS, nburn, nsave, irfperiods, optionsMS);


%% Plot
FillColor   =[.85 .85 .85];
figure

for ii = 2:size(VAR_Or,2) % Loop through VARs where target is in levels
subplot(size(Y,2),size(VAR_Or,2),ii*1)
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),[0.16])  fliplr(quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),1),'-k','LineWidth',0.5)
%     plot([quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),[0.16])' quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),[0.84])'],'-r','LineWidth',2);
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-.5 1.5])
    set(gca,'linewidth',2)
    text = [VAR_Nice{ii}];
    title(text)
    ylabel(['Labour Productivity'])
 subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),[0.16])  fliplr(quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),1),'-k','LineWidth',0.5)
    hold on
%     plot([quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),[0.16])' quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),[0.84])'],'-r','LineWidth',2);
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-0.5 1])
    set(gca,'linewidth',2)
    ylabel(['Hours'])
  subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*2+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),[0.16])...
     fliplr(quantile((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),1),'-k','LineWidth',0.5)
    hold on
%     plot([quantile((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),[0.16])' quantile((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),[0.84])'],'-r','LineWidth',2);
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylim([-2 4])
    ylabel(['Investment'])
  subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*3+ii))
  fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,4,:,1))),[0.16])...
     fliplr(quantile((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,4,:,1))),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median((squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))+squeeze(IRF_draws.(VAR_Or{ii})(:,4,:,1))),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylim([-1 2])
    ylabel(['Consumption'])
  subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*4+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,5,:,1)),[0.16])...
     fliplr(quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,5,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(squeeze(IRF_draws.(VAR_Or{ii})(:,5,:,1)),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-0.4 0.2])
    set(gca,'linewidth',2)
    ylabel(['Inflation'])
   subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*5+ii))
       fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,6,:,1)),[0.16])...
     fliplr(quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,6,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
 
    plot(median(squeeze(IRF_draws.(VAR_Or{ii})(:,6,:,1)),1),'-k','LineWidth',0.5)
    hold on
   plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    set(gca,'linewidth',2)
    ylabel(['Long Rates'])
    ylim([-0.4 0.5])
end

for ii = 1:1 % Loop through VARs where estimation is in differences for labor productivity
subplot(size(Y,2),size(VAR_Or,2),ii*1)
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))')',[0.16])  fliplr(quantile(cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))')',[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1))'),2),'-k','LineWidth',0.5)
    hold on
   plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-.5 1.5])
    set(gca,'linewidth',2)
    text = [VAR_Nice{ii}];
    title(text)
    ylabel(['Labour Productivity'])
 subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))')',[0.16])  fliplr(quantile(cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))')',[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1))'),2),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-0.5 1])
    set(gca,'linewidth',2)
    ylabel(['Hours'])
  subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*2+ii))
    fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile((cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),2)+cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),2)+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),[0.16])...
    fliplr(quantile((cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),2)+cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),2)+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median((cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),2)+cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),2)+squeeze(IRF_draws.(VAR_Or{ii})(:,3,:,1))),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-2 4])
    set(gca,'linewidth',2)
    ylabel(['Investment'])
  subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*3+ii))
   fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile((cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),2)+cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),2)+squeeze(IRF_draws.(VAR_Or{ii})(:,4,:,1))),[0.16])...
    fliplr(quantile((cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),2)+cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),2)+squeeze(IRF_draws.(VAR_Or{ii})(:,4,:,1))),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median((cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,2,:,1)),2)+cumsum(squeeze(IRF_draws.(VAR_Or{ii})(:,1,:,1)),2)+squeeze(IRF_draws.(VAR_Or{ii})(:,4,:,1))),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-1 2])
    set(gca,'linewidth',2)
    ylabel(['Consumption'])
  subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*4+ii))
   fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,5,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,5,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(squeeze(IRF_draws.(VAR_Or{ii})(:,5,:,1)),1),'-k','LineWidth',0.5)
    hold on
    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-0.4 0.2])
    set(gca,'linewidth',2)
    ylabel(['Inflation'])
   subplot(size(Y,2),size(VAR_Or,2),(size(VAR_Or,2)*5+ii))
      fill([(1:irfperiods) fliplr((1:irfperiods))],...
    [ quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,6,:,1)),[0.16])...
    fliplr(quantile(squeeze(IRF_draws.(VAR_Or{ii})(:,6,:,1)),[0.84]))],...
    FillColor,'EdgeColor','none');
    hold on
    plot(median(squeeze(IRF_draws.(VAR_Or{ii})(:,6,:,1)),1),'-k','LineWidth',0.5)
    hold on

    plot(zeros(irfperiods,1),'-k','LineWidth',0.5)
    hold off
    box on
    ylim([-0.4 0.5])
    set(gca,'linewidth',2)
    ylabel(['Long Rates'])
end

set(gcf, 'Position',  [100, 100, 1000, 800])

%% Plot FEVDs
figure 
FEVDs = [median(squeeze(FEVD_draws.Long_Run(:,1,1,:)))' median(squeeze(FEVD_draws.Max_Share(:,1,1,:)))'];
figure
plot(FEVDs(1:40,1),'-.r*','LineWidth',2)
hold on
plot(FEVDs(1:40,2),'--mo','LineWidth',2)

legend({'Long-Run','Max Share'},'fontsize',18, 'Location', 'east')
xlabel('Periods','fontsize',18)
ylabel('FEVD share','fontsize',18)
axis tight
ylim([0 1])
legend boxoff 
set(gca,'FontSize',18)
set(gcf, 'units', 'centimeters', 'Position', [1,2,25,20])


%% Extract shock decomposition - example
%HDshock(draw, variable, time, shock)

Tech_prodMS = median(squeeze(HDshock_draws.Max_Share(:, 1,: ,1))); % median contribution of tech shock to labor productivity
Init_prodMS =  median(squeeze(HDinit_draws.Max_Share(:,1,:)));
Const_prodMS =  median(squeeze(HDconst_draws.Max_Share(:,1,:)));

