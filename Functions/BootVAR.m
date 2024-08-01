function [InvA_draws,ALPHA_draws, SIGMA_draws, HDshock_draws, HDinit_draws, HDconst_draws, IRF_draws, FEVD_draws, EPS_draws] = BootVAR(X, Y, nlags, nvars, Ident, nsave, irfperiods, options)

% Returns Inverse A matrix draws, dependent on choice of identification,
% Alpha (VAR coefficients), Var-Cov matrix (Sigma), Historical shock
% decomp, initial condition contribution, constant constribution, IRFs,
% FEVD, and structural shocks (eps).


% Y and X in specification Y = Alpha*X+error. X should have constants and exogenous variables ordered last after lags of endogenous variables. 

%% Basic parameters
totiters = nsave;
T= size(Y,1);
K = size(X,2); 

%% Create priors
% First get OLS estimators
A_OLS = inv(X'*X)*(X'*Y); 
a_OLS = A_OLS(:);         

SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);   
SIGMA_OLS = SSE./(T-(K));
res = (Y - X*A_OLS);

% Storage space for posterior draws, hdecomp, FEVD, and IRFs
ALPHA_draws = nan(nsave,K,nvars);
SIGMA_draws = nan(nsave,nvars,nvars);
InvA_draws = nan(nsave, nvars, nvars);

HDshock_draws = nan(nsave,nvars,T+1,nvars);
HDinit_draws = nan(nsave,nvars,T+1);
HDconst_draws = nan(nsave,nvars,T+1);

IRF_draws = nan(nsave, nvars, irfperiods, nvars);
FEVD_draws = nan(nsave, nvars, nvars, irfperiods);
EPS_draws = nan(nsave,nvars,T);


    
    
%% Wild boostrap sampling
tic;
% disp('Number of iterations');
iters = 1;
disp("Loops completed:")
while iters < totiters  %Start the draw loop
    if mod(iters,100) == 0
        disp(iters);
        toc;
    end
    
    %Wild bootstrap: Generate new reduced form residuals with signs flipped (necessary for consistency between IV and resids, maintain positions, rather than resampling residuals across all time periods).
    
    flips = 2*(rand(T,1)>0.5)-1;
    resdraw = res.*(flips);

    if isfield(options,'IV')
        IVdraw = [flips.*options.IV(:,:)];
    end
    
    %% Initialize new data with original data
    ArtificialY = zeros(size(Y,1), size(Y,2));
    ArtificialX = zeros(size(X,1), size(X,2));
    initX = X(1,:);
    for jj = 1:size(ArtificialY,1)
        ArtificialY(jj,:) = initX*A_OLS+resdraw(jj,:);
        ArtificialX(jj,:) = initX;
        % now update X with the latest artifical draw for the first lag
        initX(end-nvars:end-1) = []; % remove last lag
        initX = [ArtificialY(jj,:) initX]; % Insert latest new Y for the next iter
    end


    % Re-estimate VAR
    ALPHA = inv(ArtificialX'*ArtificialX)*(ArtificialX'*ArtificialY);        
    SSE = (ArtificialY - ArtificialX*ALPHA)'*(ArtificialY - ArtificialX*ALPHA);   
    SIGMA = SSE./(T-(K));
    resdraw = (ArtificialY - ArtificialX*ALPHA);
    % Check stability
    Fcomp = [ALPHA(1:nvars*nlags,:)'; eye(nvars*(nlags-1)) zeros(nvars*(nlags-1),nvars)];
    if max(abs(eig(Fcomp)))>1 
       continue % Skip to next draw if fails stability test
    end
    


    
    % Choose identification and find
        PlotIRF = 0;
        switch Ident
            case 1 %Gali 1999 LR restrictions
                [InvA] = LongRun(ALPHA,SIGMA,nvars,nlags);
            case 2 % Francis et al. 2014 Max-Share identification
                [InvA] = MaxShare(ALPHA,SIGMA,nvars,nlags,options.horz, options.target); %preset to use 10 years horizon for FECD identification
            case 3 % Cholesky 
                [InvA] = chol(SIGMA,'lower');
            case 4 % Sign, zero and magnitude
                 [InvA] = LRSR_restrict_bayes(ALPHA,SIGMA,nvars,nlags, options.restrictsign,options.restrictzero, options.ratio, options.trys, options.FEVDrestrict, options.FEVDperiods);
            case 5 % Spectral identification
                [InvA] = SpecIdent(PlotIRF,ALPHA,SIGMA,nvars,nlags, options.freqlow, options.freqhigh, options.target); %Last 3 inputs (min q freq, min q freq, var input)
            case 6 % Truncated spectral identification
                [InvA] = SpecIdentLim(ALPHA,SIGMA,nvars,nlags, options.freqlow, options.freqhigh, options.target, options.horzlim);
            case 7 % Restrict impact matrix
                [InvA] = A0restrict(ALPHA,SIGMA,nvars,nlags, options.trys);
            case 8 % IV var
                [InvA] = IVrestrict(ALPHA,SIGMA,nvars,nlags, IVdraw, resdraw);
            otherwise
                disp('You have not picked an identification methodology!')
        end
        
        % Rotate a shock to be positive or negative if needed
        if isfield(options,'rotatesign')==1
            if options.rotatesign==1
                if InvA(1,1)<0
                    InvA(:,1) = -InvA(:,1);
                end
            end
        end
        
        % Save draws
        if ~isnan(InvA)
            InvA_draws(iters,:,:)  = InvA;
            ALPHA_draws(iters,:,:) = ALPHA;
            SIGMA_draws(iters,:,:) = SIGMA;
            
            % Make structural shocks series           
            EPS_draws(iters,:,:) = InvA\res';  
            warning('off','MATLAB:singularMatrix')
            %HDshock(variable, time, shock) contributions - %HDinit and const are (var, time).
            [HDshock_draws(iters,:,:,:), HDinit_draws(iters,:,:), HDconst_draws(iters,:,:)]...
            = HDecomp(InvA,ALPHA,nvars,nlags,T, Y, X);
                
            %Get IRFs
            [IRF_draws(iters,:,:,:)] = IRFrun(InvA,ALPHA,nvars,nlags,irfperiods);
            
            %Get FEVD - calculate contribution for each variable - if
            %options.diff = 1 is chosed, the function will cumulate the
            %FEVD as if the estimation was on the variable in differences,
            %but the contribution to the variance of the variable in levels
            %is desired.
            % FEVD(endogvar,shock ,time)
            if isfield(options, 'diff') && options.diff == 1
                [FEVD_draws(iters,:,:,:)] = FEVDdiff(InvA,ALPHA,nvars,nlags,irfperiods,options.target); %Horizon of FEVD
            else
                [FEVD_draws(iters,:,:,:)] = FEVD(InvA,ALPHA,nvars,nlags,irfperiods);
            end
           
        end
        
        iters = iters+1;

      
end

end
