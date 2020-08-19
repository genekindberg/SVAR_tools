function [InvA_draws,ALPHA_draws, SIGMA_draws, HDshock_draws, HDinit_draws, HDconst_draws, IRF_draws, FEVD_draws, EPS_draws] = BVAR_Minn(X, Y, nlags, nvars, Ident, nburn, nsave, irfperiods, options)

totiters = nburn+nsave;
T= size(Y,1);
K = size(X,2); % 1 for constant
%Based on Koop and Korobolis code
%% Create priors
% First get ML estimators
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of parameters, i.e. it holds
                          % that a_OLS = vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);   % Sum of squared errors
SIGMA_OLS = SSE./(T-(K));

% Initialize Bayesian posterior parameters using OLS values
alpha = a_OLS;     % This is the single draw from the posterior of alpha
ALPHA = A_OLS;     % This is the single draw from the posterior of ALPHA Just reshape previous alpha vector
SSE_Gibbs = SSE;   % This is the single draw from the posterior of SSE
SIGMA = SIGMA_OLS; % This is the single draw from the posterior of SIGMA

% Storage space for posterior draws and hdecomp
ALPHA_draws = nan(nsave,K,nvars);
SIGMA_draws = nan(nsave,nvars,nvars);
InvA_draws = nan(nsave, nvars, nvars);

HDshock_draws = nan(nsave,nvars,T+1,nvars);
HDinit_draws = nan(nsave,nvars,T+1);
HDconst_draws = nan(nsave,nvars,T+1);

IRF_draws = nan(nsave, nvars, irfperiods, nvars);
FEVD_draws = nan(nsave, nvars, nvars, irfperiods);
EPS_draws = nan(nsave,nvars,T);

%% hyper -params for minnesota prior

A_prior = [0.9*eye(nvars); zeros((nlags-1)*nvars,nvars);zeros(size(X,2)-nlags*nvars,nvars)];  %<---- prior mean of ALPHA (parameter matrix) AR1 is 0.9 and others are 0.
a_prior = A_prior(:); %<---- prior mean of alpha (parameter vector)
    
%tightness priors
a_bar_1 = 1;%0.5; %tightness on own lags
a_bar_2 = 1;%1; %tightness on lags of other variables
a_bar_3 = 1; %1 % tightness on own lags - how much tighter prior gets over time
a_bar_4 = 10^2; %tightness on constant

sigma_sq = zeros(nvars,1); % vector to store residual variances
%Set prior for VAR variance as variance of errors from AR equation - get
%variances from ARs first
for i = 1:nvars
        % Create lags of dependent variable in i-th equation
        Ylag_i = mlag2(Y(:,i),nlags);
        Ylag_i = Ylag_i(nlags+1:T,:);
        % Dependent variable in i-th equation
        Y_i = Y(nlags+1:T,i);
        % OLS estimates of i-th equation
        alpha_i = inv(Ylag_i'*Ylag_i)*(Ylag_i'*Y_i);
        sigma_sq(i,1) = (1./(T-nlags+1))*(Y_i - Ylag_i*alpha_i)'*(Y_i - Ylag_i*alpha_i);
end

    % Variance on priors
    V_i = zeros(K,nvars);
    
    % index in each equation which are the own lags
%     for ii=1:n
%         for jj=1:p
%             phi0((jj-1)*n+ii,(jj-1)*n+ii)=(1/arvar(ii,1))*(lambda1/jj^lambda3)^2;
%         end
%     end
    ind = zeros(nvars,nlags);
    for i=1:nvars
        ind(i,:) = i:nvars:nvars*nlags; %constant is at the end!
    end
    for i = 1:nvars  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            lagord = ceil(j/nvars);
                if j>nlags*nvars
                    V_i(j,i) = a_bar_4*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = (1/sigma_sq(i,1))*(a_bar_1./(lagord^a_bar_3))^2; % variance on own lags           
                else
                    for kj=1:nvars
                        if find(j==ind(kj,:))>0 % Find which variable the current j belongs to
                            ll = kj;                   
                        end
                    end                 % variance on other lags   
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((nlags^2)*sigma_sq(ll,1));      
                end
        end
    end
    
    % Now V is a diagonal matrix with diagonal elements the V_i
    V_prior = diag(V_i(:));  % this is the prior variance of the vector alpha
    
    % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior)) -Var(1)
    v_prior = nvars;
    S_prior = SIGMA_OLS;
    inv_S_prior = inv(S_prior); 
    
    
%% Gibbs sampler
tic;
% disp('Number of iterations');
for iters = 1:totiters  %Start the Gibbs "loop"
%     if mod(iters,200) == 0
%         disp(iters);
%         toc;
%     end
    
    %Draw posterior beta, V and sigma
    V_post = inv(inv(V_prior) + kron(inv(SIGMA),X'*X));
    a_post = V_post*(inv(V_prior)*a_prior + kron(inv(SIGMA),X'*X)*a_OLS);
    alpha = a_post + chol(nspd(V_post))'*randn(K*nvars,1); % Draw of alpha
    
    ALPHA = reshape(alpha,K,nvars); % Draw of ALPHA
    %% Check stbility and draw again if needed
    Fcomp = [ALPHA(1:nvars*nlags,:)'; eye(nvars*(nlags-1)) zeros(nvars*(nlags-1),nvars)];
    stopcount = 1;
    while max(abs(eig(Fcomp)))>1
        alpha = a_post + chol(nspd(V_post))'*randn(K*nvars,1);
        ALPHA = reshape(alpha,K,nvars);
        Fcomp = [ALPHA(1:nvars*nlags,:)'; eye(nvars*(nlags-1)) zeros(nvars*(nlags-1),nvars)];
        stopcount=stopcount+1;
        if stopcount>1000
            break
        end
    end
    
    % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
    v_post = T + v_prior;
    S_post = S_prior + (Y - X*ALPHA)'*(Y - X*ALPHA);
    A = chol(nspd(inv(S_post)))'*randn(size(inv(S_post),1),v_post);
    A = A*A';
    SIGMA = nspd(inv(A));% Draw SIGMA - stabilize to be pos def
    if iters > nburn   

    
    % Choose identification and find
        PlotIRF = 0;
        switch Ident
             case 1 %Gali 1999 LR restrictions
                [InvA] = gali1999(PlotIRF,ALPHA,SIGMA,nvars,nlags);
             case 2 % Francis max FEVD ident
                [InvA] = Francis2014(PlotIRF,ALPHA,SIGMA,nvars,nlags,options.horz, options.target); %preset to use 10 years horizon for FECD identification
             case 3 % Fisher IST and neutral shock ID
                [InvA] = Fisher2006_GMM(PlotIRF,ALPHA,SIGMA,nvars,nlags);
             case 4 % Identification max variance at low freq
                [InvA] = SpecIdent(PlotIRF,ALPHA,SIGMA,nvars,nlags, options.freqlow, options.freqhigh, options.target); %Last 3 inputs (min q freq, min q freq, var input)
            case 5 % SR and LR restriction identificaiton
                [InvA] = LRSR_restrict_bayes(ALPHA,SIGMA,nvars,nlags, options.restrictsign,options.restrictzero, options.restrict_horiz, options.trys);
            case 6 % Modified max share
                [InvA] = ModMaxShare(PlotIRF,ALPHA,SIGMA,nvars,nlags,options.horz, options.target);
            case 7 % Limited spectral
                [InvA] = SpecIdentLim(PlotIRF, ALPHA,SIGMA,nvars,nlags, options.freqlow, options.freqhigh, options.target, options.horzlim);
            case 8 % Cholesky 
                [InvA] = chol(SIGMA,'lower');
            otherwise
                disp('You have not picked an identification methodology!')
        end
        
        if isfield(options,'rotatesign')==1
            if options.rotatesign==1
                if InvA(1,1)<0
                    InvA(:,1) = -InvA(:,1);
                end
            end
        end
        
        if ~isnan(InvA)
            InvA_draws(iters-nburn,:,:)  = InvA;
            ALPHA_draws(iters-nburn,:,:) = ALPHA;
            SIGMA_draws(iters-nburn,:,:) = SIGMA;
            % Make structural shocks series
            
             res = (Y - X*ALPHA);
            EPS_draws(iters-nburn,:,:) = InvA\res';  
%             
%  

 
        
                %Shock decomp
                %HDshock(variable, time, shock) contributions - %HDinit and const =
                %1 by time.
            [HDshock_draws(iters-nburn,:,:,:), HDinit_draws(iters-nburn,:,:), HDconst_draws(iters-nburn,:,:)]...
            = HDecomp(InvA,ALPHA,nvars,nlags,T, Y, X);
                
            %Get IRFs
            [IRF_draws(iters-nburn,:,:,:)] = IRFrun(InvA,ALPHA,nvars,nlags,irfperiods);
            %Get FEVD
            if options.diff == 1
                [FEVD_draws(iters-nburn,:,:,:)] = FEVDdiff(InvA,ALPHA,nvars,nlags,irfperiods,options.target); %Horizon of FEVD
            else
                [FEVD_draws(iters-nburn,:,:,:)] = FEVD(InvA,ALPHA,nvars,nlags,irfperiods);
            end
           
        end
        
       
    end
      
end

end
