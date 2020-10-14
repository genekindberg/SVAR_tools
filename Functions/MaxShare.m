function [InvA]=MaxShare(B,varmat,nvars,nlags,horzlim, target)

M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

%Take cholesky decomp of var-cov error mat
Atilde = chol(varmat)';

%Selection matrix for the variable targetted for max-share
ei = zeros(nvars,1); 
ei(target,1) = 1;      

% calculate impact mat at each horizon up to horzlim
Bi = zeros(1,nvars,horzlim+1);
for l=0:horzlim
    C=M^l;
    Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars); 
end

%compute V - FEV
V = zeros(nvars,nvars);

for l=0:horzlim 
    V = V + (Bi(:,:,l+1)*Atilde)'*(Bi(:,:,l+1)*Atilde);
end

VHat = V(1:nvars,1:nvars);

[eigenVector,eigenValue]=eig(VHat); 
%sort so that largest eigenvector is at the top
sorted=[diag(eigenValue) (1:1:(nvars))'];
order=sortrows(sorted,-1);
VHatsort = eigenVector(:,order(:,2)');

InvA = Atilde*VHatsort; %Create new ident matrix


end
