function [IRFs] = IRF_run(InvA,ALPHA,nvars,nlags,periods)

% Create companion  matrix
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=ALPHA(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

IRFs = zeros(nvars, periods,nvars);
for ii = 1:nvars
gamma = zeros(nvars,1);
gamma(ii,1) = 1;
alpha=InvA*gamma;
Impulse=[alpha; zeros(nvars*nlags-nvars,1)];

    for k=1:periods
    %variable, time, shock
        temp = (M^(k-1)*Impulse);
    %variable time shock
        IRFs(:,k,ii)=temp(1:nvars,:);
    end
end

end
