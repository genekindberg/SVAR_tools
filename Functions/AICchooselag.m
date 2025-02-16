function optimalLag = AICchooselag(Data,maxlag)
    [nObs, nvars] = size(Data); % Number of observations and variables
    aicValues = zeros(maxlag, 1); % Initialize AIC values
    for nlags = 1:maxlag
        % Prepare the lagged data
            X = [];
            for p=1:nlags
                X(:,1+(p-1)*nvars:p*nvars)=Data((nlags+1-p):(nObs-p),:);
                %gives (T-nlags) x nvars*nlags matrix of lags for all variables                                         
            end
            % Add constant
            const = ones(nObs-nlags,1);

            X = [X const];
            %rescaling Y variables since we loose nlags observations
            Y=Data((nlags+1):end,:);

        
            % OLS estimation: b = (X'X)^(-1)X'Y
            b = (X' * X) \ (X' * Y); % Coefficients estimation
            
            % Calculate residuals
            Y_hat = X * b; % Fitted values
            residuals = Y - Y_hat; % Residuals
            
            % Calculate log-likelihood
            logLikelihood = -0.5 * nObs * nvars * log(2 * pi) - 0.5 * nObs * log(det(residuals' * residuals / nObs));
            
            % Number of parameters: k * (p * n + 1)
            numParams = nvars * (nlags * nvars + 1);
        
            % Calculate AIC
            aicValues(nlags) = -2 * logLikelihood + 2 * numParams;
    end

    % Find the optimal lag length (the one with the minimum AIC value)
    [~, optimalLag] = min(aicValues);

    X = sprintf('Optimal number of lags (AIC): %d',optimalLag);
    disp(X)

end