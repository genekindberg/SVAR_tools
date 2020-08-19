function [X, Y] = MakeData(outputinit, nlags)

            data = outputinit;
            X = [];
           [T,nvars]=size(data);
                for p=1:nlags
                    X(:,1+(p-1)*nvars:p*nvars)=data((nlags+1-p):(T-p),:);
                    %gives (T-nlags) x nvars*nlags matrix of lags for all variables                                         
                end

            const = ones(T-nlags,1);
            X = [X const];                
            %rescaling variables since we lose nlags observations through the lagging
            Y=data((nlags+1):T,:);
            
end