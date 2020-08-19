function [InvA] = gali1999(PlotIRF, B,varmat,nvars,nlags)



% Create companion  matrix
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

% Create inverse A matrix based on long run restrictions - as in Blachard
% Quah 1989.
Finf = inv(eye(length(M))-M); % from the companion
Finf = Finf(1:nvars,1:nvars);
D  = chol(Finf*varmat*Finf')'; % identification: u2 has no effect on y1 in the long run
InvA = Finf\D;

if PlotIRF == 1
%variable time shock
periods = 12;
IRF = IRFrun(InvA,B,nvars,nlags,periods);

figure
subplot(2,2,1)
plot(1:periods,cumsum(IRF(1,:,1)),'-','LineWidth',2)

xlabel('Productivity')
title('Technology shock')
axis tight
ylim([0.4 1.1])
subplot(2,2,2)
summed = cumsum(IRF(1,1:end,2));
plot(1:periods,[summed(1,2:end),0.01],'-','LineWidth',2)
hold on
plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
xlabel('Productivity')
title('Non-technology shock')
axis tight
hold off
ylim([-0.2 0.6])
subplot(2,2,3)
plot(1:periods,cumsum(IRF(2,:,1)),'-','LineWidth',2)
hold on
plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
hold off
xlabel('Hours')
title('Technology shock')
axis tight
ylim([-1 0.6])
subplot(2,2,4)
plot(1:periods,cumsum(IRF(2,:,2)),'-','LineWidth',2)
hold on
plot(1:periods,[zeros(1,periods)],'-k','LineWidth',0.5)
hold off
xlabel('Hours')
title('Non-technology shock')
axis tight
ylim([0.4 2])
set(gcf, 'units', 'centimeters', 'Position', [1,2,20,12])

end


end


