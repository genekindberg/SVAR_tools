function [HDshock, HDinit, HDconst] = HDecomp(InvA,B,nvars,nlags,nobs, Y, X)

% outputs:         %HDshock(variable, time, shock) contributions - %HDinit and const = [1 ,time].
res = (Y - X*B);

% Make companion matrix
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

% Make structural shocks series
eps = InvA\res';

% Make structural errors series
invA_big = zeros(nvars*nlags,nvars);
invA_big(1:nvars,:) = InvA;


%% Shock Decomp
Icomp = [eye(nvars) zeros(nvars,(nlags-1)*nvars)];
HDshock_big = zeros(nlags*nvars,nobs+1,nvars); % used in calculation (requires previous values in matrix mult.)
HDshock = zeros(nvars,nobs+1,nvars);
    for j=1:nvars % for each variable
        eps_big = zeros(nvars,nobs+1); % matrix of shocks conformable with companion
        eps_big(j,2:end) = eps(j,:);
        for i = 2:nobs+1
            HDshock_big(:,i,j) = invA_big*eps_big(:,i) + M*HDshock_big(:,i-1,j);
            HDshock(:,i,j) =  Icomp*HDshock_big(:,i,j); % remove lagged values
        end
    end
    
%% Initial value
HDinit_big = zeros(nlags*nvars,nobs+1);
HDinit = zeros(nvars, nobs+1);
HDinit_big(:,1) = X(1,1:nlags*nvars)';
HDinit(:,1) = Icomp*HDinit_big(:,1);
    for i = 2:nobs+1
        HDinit_big(:,i) = M*HDinit_big(:,i-1);
        HDinit(:,i) = Icomp *HDinit_big(:,i);
    end
    
% Constant
HDconst_big = zeros(nlags*nvars,nobs+1);
HDconst = zeros(nvars, nobs+1);
CC = zeros(nlags*nvars,1);
DD = zeros(nlags*nvars,1);
CC(1:nvars,:) = B(nvars*nlags+1,:)';
if size(X,2)<nvars*nlags+2
    for i = 2:nobs+1
         HDconst_big(:,i) = CC + M*HDconst_big(:,i-1);
         HDconst(:,i) = Icomp * HDconst_big(:,i);
    end
else
    for i = 2:nobs+1
         DD(1:nvars,:) = (X(i-1,nlags*nvars+2:end).*B(nvars*nlags+2:end,:))';
         HDconst_big(:,i) = CC +DD+ M*HDconst_big(:,i-1);
         HDconst(:,i) = Icomp * HDconst_big(:,i);
    end
end
    
    
sumshocks = squeeze(sum(HDshock,3));
alldecomp =(sumshocks(:,2:end)'+HDinit(:,2:end)'+HDconst(:,2:end)');
if abs(Y-(sumshocks(:,2:end)'+HDinit(:,2:end)'+HDconst(:,2:end)')) >.0001;
    disp('Shocks do not sum to variable')
end


end
