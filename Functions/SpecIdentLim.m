function [InvA] = SpecIdentLim(B,varmat,nvars,nlags, a, b, varmax, horzlim)

M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

%Take cholesky decomp of 
Atilde = chol(varmat)';


N=200; % Number divisions over 2*pi span to collect for Reimann sum
kk=1:N/2; % Do riemann sum over N components of 2*pi 
%% Calculate variance at different frequencies
ei = zeros(nvars,1); 
ei(varmax,1) = 1;  % Selection matrix for variable of interest
%Get variance at desired frequencies frequencies
aa=floor((2*pi/a)/(pi/(N/2))); 
bb=max(floor((2*pi/b)/(pi/(N/2))),1); %so does not go to zero at very low frequencies
count = 1;
for ll = bb:aa
    
    k=kk(ll);
    omeg=2*pi*k/N;
    i=sqrt(-1);
    z=exp(-i*omeg);
      
    % Make impact matrix at each horizon
    Bi = zeros(1,nvars,horzlim+1);
    for l=0:horzlim
        C=M^l;
        Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars); % Bi are VMA coeffs on how ith variable responds to different shocks l+1 periods back                                     % Bi(:,:,l+1) is the coef on epsilon(t-l)
    end
    
    % now make cov Y(t)*Y(t-j) for j=0:horzlim - assuming that each y can
    % be represented by its ma(horzlim) representation, as opposed to the
    % infinity MA representation y = inv(I-A)eps
    % each period j covariance can be represented as Y(t)Y(t-10) = Bi0*V*Bi10
    % Y(t)Y(t-9) = Bi0*V*Bi9+ Bi1*V*Bi10 etc....
    % Y(t)Y(t+9) = Bi9*V*Bi0+ Bi10*V*Bi1
    %% for Y(t)Y(t-ii)
    A=Bi(:,:,1);
    if nlags > 1
        for ii = 2:horzlim+1
            A=A+Bi(:,:,ii)*(z^ii);
        end
    end
    
    
    Sy(:,:,count)=(A*Atilde)'*(A*Atilde);
    count = count+1;

end

%desired frequency
%freq= [0:pi/(N/2):pi]';
%frequency of lower and upper band 2*pi/a 2*pi/b*(pi/(N/2))

Spect=zeros(size(Sy(:,:,1),1));
%Summable frequencies
for jj = 1:size(Sy,3)
    Spect = Spect+real(Sy(:,:,jj));
end

%% find Max explained variance shock
%Get max explained using shock at this frequency
% Select all shocks - do not exclude first shock
VHat = Spect;

[eigenVector,eigenValue]=eig(VHat); 
%sort so that largest eigenvector is at the top
sorted=[diag(eigenValue) (1:1:(nvars))'];
order=sortrows(sorted,-1);
VHatsort = eigenVector(:,order(:,2)');

InvA = Atilde*VHatsort;
%Flip signs on first shock if negative
% Check sign at period 40 to check sign and whether flip
Impulse = zeros(nvars,1); Impulse(1,1)=1;
Impulse = InvA*Impulse;
Impulse = [Impulse; zeros(nvars*nlags-nvars,1)];
temp = (M^(a-1)*Impulse);
%variable time shock


if sum(temp(varmax,1))<0
    InvA(:,1) = -InvA(:,1); % Flip sign so that impact on prod is positive
end

% if InvA(varmax,2)<0
%     InvA(:,2) =  -InvA(:,2);
% end




end


