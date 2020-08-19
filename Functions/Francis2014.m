function [InvA]=Francis2014(PlotIRF, B,varmat,nvars,nlags,horzlim, target)

%make companion-form matrix.  Excludes constant term (Last in the 'B' matrix).
%Y(t;t-1;t-2) = M*Y(t-1;t-2;t-3)
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
    Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars); % Bi are VMA coeffs on how ith variable responds to different shocks l+1 periods back
                                       % Bi(:,:,l+1) is the coef on epsilon(t-l)
end

%compute V - FEV
V = zeros(nvars,nvars);

%impact of all shocks on Adj TFP
for l=0:horzlim 
    V = V + (Bi(:,:,l+1)*Atilde)'*(Bi(:,:,l+1)*Atilde);
end

% Select all shocks - do not exclude first shock as in BS2011
VHat = V(1:nvars,1:nvars);

[eigenVector,eigenValue]=eig(VHat); 
%sort so that largest eigenvector is at the top
sorted=[diag(eigenValue) (1:1:(nvars))'];
order=sortrows(sorted,-1);
VHatsort = eigenVector(:,order(:,2)');

InvA = Atilde*VHatsort; %Create new ident matrix

if InvA(target,1)<0
    InvA(:,1) = -InvA(:,1); % Flip sign so that impact on prod is positive
end

% Plot IRFs as in Francis et al. 2011
if PlotIRF == 1
%variable time shock
periods = 25;
IRF = IRFrun(InvA,B,nvars,nlags,periods);
%prod hours, hours per cap, GDP, Consumption/GDP,invest/GDP

figure
subplot(5,1,1)
%note - IRF=(endog var, time, shock) ordering
plot(1:periods,IRF(1,:,1),'-','LineWidth',2)
ylabel('Productivity')
title('Technology shock')
axis tight
 ylim([0 0.8])
subplot(5,1,2)
plot(1:periods,[IRF(2,:,1)],'-','LineWidth',2)
% % hold on
% % plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% % hold off
ylabel('Hours')
axis tight
 ylim([-0.5 1])

subplot(5,1,3)
plot(1:periods,[IRF(1,:,1)+IRF(2,:,1)],'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('GDP')
axis tight
 ylim([0 1])

subplot(5,1,4)
plot(1:periods,IRF(3,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Consumption')
axis tight
 ylim([-0.5 1])

subplot(5,1,5)
plot(1:periods,IRF(4,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Investment')
axis tight
ylim([-1 3.5])
set(gcf, 'units', 'centimeters', 'Position', [4,-10,12,25])
end


end
