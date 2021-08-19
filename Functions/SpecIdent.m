function [InvA] = SpecIdent(PlotIRF, B,varmat,nvars,nlags, a, b, varmax)

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
%Get variance at different frequencies
aa=floor((2*pi/a)/(pi/(N/2))); 
bb=max(floor((2*pi/b)/(pi/(N/2))),1); %so does not go to zero at very low frequencies
count = 1;
for ll = bb:aa
    
    k=kk(ll);
    omeg=2*pi*k/N;
    i=sqrt(-1);
    z=exp(-i*omeg);
      
    % To make (I-A1e^-iw - A2e^-iw2....)^(-1)
    A=B(1:nvars,:)'*z;
    if nlags > 1
        for ii = 2:nlags
            A=A+B(nvars*(ii-1)+1:nvars*(ii),:)'*(z^ii);
        end
    end
    DD = ei'*inv(eye(nvars)-A); %LR impact on variable I (using selection mat ei)
    
    
    %Spectrum of first variable 
%     Sy(:,:,ll)=(DD*Atilde)'*(Atilde'*ctranspose(DD))';
    Sy(:,:,count)=(DD*Atilde)'*(DD*Atilde);
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


if temp(varmax,1)<0
    InvA(:,1) = -InvA(:,1); % Flip sign so that impact on prod is positive
end
%%  calculate IRF
if PlotIRF == 1
%variable time shock
periods = 20;
IRF = IRFrun(InvA,B,nvars,nlags,periods);
% productivity, hours, GDP, Consumption, Investment

figure
subplot(3,4,1)
plot(1:periods,IRF(1,:,1),'-','LineWidth',2)
ylabel('GDP')
title('Output')
axis tight
ylim([-0.5 0.9])
subplot(3,4,2)
plot(1:periods,[IRF(2,:,1)],'-','LineWidth',2)
% % hold on
% % plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% % hold off
ylabel('Hours')
axis tight
ylim([-0.75 0.75])

subplot(3,4,3)
plot(1:periods,IRF(3,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Investment')
axis tight
ylim([-1 3])

subplot(3,4,4)
plot(1:periods,IRF(4,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Consumption')
axis tight
ylim([-0.5 1])

subplot(3,4,5)
plot(1:periods,IRF(5,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Unemployment')
axis tight
ylim([-.75 .75])

subplot(3,4,6)
plot(1:periods,IRF(6,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Inflation')
axis tight
ylim([-.75 .75])

subplot(3,4,7)
plot(1:periods,IRF(7,:,1),'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('FedFunds')
axis tight
ylim([-.75 .75])

subplot(3,4,8)
plot(1:periods,[IRF(7,:,1)-IRF(6,:,1)],'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Real interest')
axis tight
ylim([-.75 .75])

subplot(3,4,9)
plot(1:periods,[IRF(8,:,1)],'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Fernald TFP')
axis tight
ylim([-.75 .75])

subplot(3,4,10)
plot(1:periods,[IRF(1,:,1)-IRF(2,:,1)],'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('L Prod')
axis tight
ylim([-.75 .75])

subplot(3,4,11)
plot(1:periods,[IRF(10,:,1)],'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Real Wage')
axis tight
ylim([-.75 .75])

subplot(3,4,12)
plot(1:periods,[IRF(11,:,1)],'-','LineWidth',2)
% hold on
% plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
% hold off
ylabel('Labour Share')
axis tight
ylim([-.75 .75])

set(gcf, 'units', 'centimeters', 'Position', [4,-10,35,20])

end



end


