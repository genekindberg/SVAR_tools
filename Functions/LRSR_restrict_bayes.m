function [InvA] = LRSR_restrict_bayes(B,varmat,nvars,nlags, restrictsignI,restrictzeroI, trys, FEVDrestrict,FEVDperiods)

if nargin < 9
    FEVDrestrict=[];
    FEVDperiods=[];
end
% B = coefficient mat (const in final column)
% varmat        % cov mat of residuals
% nvars, nlags  % number of endog vars, lags
% restrictsignI  % matrix containing separate array for each shock impact and
% horizon: % [Shock, Variable, Period, Pos(1)/Neg(-1), Magnitude (at least x response)]
% restrictzeroI: Matrix containing
% inf reflects infinite horizon - '1' is impact matrix.
% trys          % number of trys to match sign restriction - gives up and
% returns nan matrix if loops exceed this.


%% Store of successful draws

InvA = nan(nvars,nvars); % ensures that at least a NaN matrix is returned if fails to detect matching shock matrix.

%% Create companion mat and calculate relevant IRFs
if restrictzeroI ~isempty(restrictzeroI);
    all_horz = unique([restrictsignI(:,3);restrictzeroI(:,3)]); % Find all unique horizons to investigate for sign and zero restrictions
else 
    all_horz = unique([restrictsignI(:,3)]);
end
% Create companion coefficient matrix.
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

%Take cholesky decomp of varmat
Atilde = chol(varmat)';

%Calculate IRFs for each required horizon and stack
inc_inf = find(all_horz==inf);
if inc_inf>0
    all_horz(inc_inf) = [];
end

IRF = []; % for stacking IRFs for each horizon
% Create 'L' impact matrix for each horizon restriction
for k=1:size(all_horz,2)
    %variable, time, shock
        temp = (M^(all_horz(k)));
    %variable time shock
        IRF=[IRF ; temp(1:nvars,1:nvars)*Atilde];
end

%add inf horizon impact to end
if inc_inf>0
    Finf = inv(eye(length(M))-M); % from the companion
    Finf = Finf(1:nvars,1:nvars)*Atilde;
    IRF = [IRF; Finf];
    all_horz = [all_horz;inf];
end


% Create restriction and magnitude matrices
restrictsign = cell(1,nvars);
magnitudes = cell(1,nvars);
restrictzero = cell(1,nvars);

orderlookup = [all_horz (1:size(all_horz))' ];
% Set up sign restriction matrices
for ii = 1:nvars
    temp = restrictsignI(restrictsignI(:,1)==ii,:);
    if ~isempty(temp)
        for jj = 1:size(temp,1)
            initial = zeros(1,nvars*size(all_horz,1));
            time = orderlookup(find(orderlookup==temp(jj,3)),2); % Find what IRF set is member of
            initial(1,(time-1)*nvars+temp(jj,2))=temp(jj,4); % Set equal to positive or negative 
            restrictsign{ii} = [restrictsign{ii};initial];
            magnitudes{ii} = [magnitudes{ii};temp(jj,5)]; % Magnitude size
        end
    end
end

%Set up zero restriction matrices
if restrictzeroI ~isempty(restrictzeroI);
    for ii = 1:nvars
        temp = restrictzeroI(restrictzeroI(:,1)==ii,:);
        if ~isempty(temp)
            for jj = 1:size(temp,1)
                initial = zeros(1,nvars*size(all_horz,1));
                time = orderlookup(find(orderlookup==temp(jj,3)),2); % Find what IRF time period set is member of
                initial(1,(time-1)*nvars+temp(jj,2))=1; % Set equal to positive or negative 
                restrictzero{ii} = [restrictzero{ii};initial];
            end
        end
    end
else
    for ii = 1:nvars
        restrictzero{ii} = [];
    end
end
        
            
        


for iter = 1:trys

%% Impose zero restriction
%sort shocks by those with most restrictions first - otherwise may not be
%able to find nullspace of some shocks
% if there are no zero restrictions, will simply draw an orthonormal
% matrix.
    order = [];
    for bb= 1:size(restrictzero,2)
        tem = size(restrictzero{bb},1);
        order = [order tem];
    end
    order = [order' (1:nvars)'];
    order= sortrows(order,-1);

    Q = []; % holding candidate draw imposing zero restrictions
    for ii=order(:,2)'
        if ~isempty(restrictzero{ii})
            R=[restrictzero{ii}*IRF;Q'];
            N=null(R)';
            X = randn(nvars,1);
            q=N'*(N*X./norm(N*X));
            Q = [Q q];
        else
            if isempty(Q)
                X = randn(nvars,1);
                q = X./norm(X);
                Q = [Q q];
            else

                R=[Q'];
                N=null(R)';
                X = randn(nvars,1);
                q=N'*(N*X./norm(N*X));
                Q = [Q q];
            end
        end
    end
    %Reorder Q according to desired shock order
    reorder = [];
    for jj = 1:nvars
        reorder = [reorder find(jj==order(:,2)')];
    end

    Q = Q(:,reorder);
    %% Check sign restrictions hold on zero rest mat.- add new matrix if successful


    success = 1;

    % Check the sign restrictions on Q conditional on zero restrictions - stop loop early if individual restriction invalid:
    for jj = 1:nvars
        if isempty(restrictsign{jj})
            continue
        end
        if restrictsign{jj}*IRF*Q(:,jj)>magnitudes{jj}
            continue
        elseif restrictsign{jj}*IRF*-Q(:,jj)>magnitudes{jj}
            Q(:,jj) = -Q(:,jj);
            continue
        else
            success = 0;
            break
        end
    end
    %Go to next loop if sign restrict failure
    if success == 0
        continue
    end

    %Check zero restrictions again if sign restrictions confirmed
    for jj = 1:nvars
        if isempty(restrictzero{jj})
            continue
        end
        if abs(restrictzero{jj}*IRF*Q(:,jj))<1e-5
            continue
        else
            success = 0;
            break
        end
    end
    
    %Check FEVD restrictions if all other restrictions confirmed
    if isempty(FEVDrestrict)
        ;
    else
        %FEVD variable is written as (endogvar,shock ,time)
        %FEVDrestrict is written as {time period}(endogenous variable, shock with largest
        %impact, shock with smaller imapact)
        [FEVDcheck] = FEVD(Atilde*Q,B,nvars,nlags,max(FEVDperiods));
        for ii=1:size(FEVDperiods,2)
            for jj=1:size(FEVDrestrict,1)
                if FEVDcheck(FEVDrestrict{ii}(jj,1),FEVDrestrict{ii}(jj,2),FEVDperiods(ii))>FEVDcheck(FEVDrestrict{ii}(jj,1),FEVDrestrict{ii}(jj,3),FEVDperiods(ii))
                    continue
                else
                    success=0;
                end
            end
        end
    end
    
    if success==1
        InvA=Atilde*Q;
        break
    end


end


