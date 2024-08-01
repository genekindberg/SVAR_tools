function [InvA]=IVrestrict(B,varmat,nvars,nlags, IV, eps)

% reduced form residuals to be instrumented:
epsP = eps(:,1);
epsQ = eps(:,2:end); % residuals for non-instrumented variables
% Get location of first and last non-nan to make comparable IV reg
FirstL = find(~isnan(IV), 1, 'first');
LastL = find(~isnan(IV), 1, 'last');
% Common sample for IV and instrumented variables shock
epsPsamp = epsP(FirstL:LastL,1);
epsQsamp = epsQ(FirstL:LastL,:);
IVsamp = IV(FirstL:LastL,1);
%Reg resid on IV and get fitted values
BetaIV = (IVsamp'*IVsamp)\(IVsamp'*epsPsamp);
epsFit = BetaIV*IVsamp;

InvAiv = nan(nvars,nvars);
% Now do second stage of 2sls - regress non-IV RF resids on fitted IV resid
sq_spratio = zeros(size(epsQ,2),1);
for ii=1:nvars-1
  Betatemp = (epsFit'*epsFit)\(epsFit'*epsQsamp(:,ii));
  sq_spratio(ii) = Betatemp;
end

%follow equations in footnote 4 of Gertler Karadi 2015
% calculate new error matrix
allerrors = [epsPsamp epsQsamp];
Sigma = (1/(size(allerrors,1)-nvars*nlags-nvars))*allerrors'*allerrors;
Sig_11 = Sigma(1,1);
Sig_21 = Sigma(2:end,1);
Sig_22 = Sigma(2:end,2:end);

Q = sq_spratio*Sig_11*sq_spratio'-(Sig_21*sq_spratio'+sq_spratio*Sig_21')+Sig_22;
s12s12 = (Sig_21-sq_spratio*Sig_11)'*inv(Q)*(Sig_21-sq_spratio*Sig_11);
sp  = sqrt(Sig_11 - s12s12);

InvA = zeros(nvars, nvars);
InvA(:,1) = [1; sq_spratio]*sp;


end
