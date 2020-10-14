function [InvA] = LongRun(B,varmat,nvars,nlags)



% Create companion  matrix
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

% Create inverse A matrix based on long run restrictions - as in Blanchard
% Quah 1989.
Finf = inv(eye(length(M))-M); % from the companion
Finf = Finf(1:nvars,1:nvars);
D  = chol(Finf*varmat*Finf')'; 
InvA = Finf\D;


end


