function [FEVDcont] = FEVDdiff(InvA,B,nvars,nlags,horzlim, target)
%FEVD(:,:,:) = (endogvar,shock ,time)
%make companion-form matrix.  Excludes constant term (at end of beta matrix).
%Y(t;t-1;t-2) = M*Y(t-1;t-2;t-3)
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=B(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);
%Set struct shock matrix
Atilde = InvA;
FEVD = nan(nvars,nvars,horzlim);
for ii=1:nvars     
% calculate coefficient lag polynomial
    ei = zeros(nvars,1); 
    ei(ii,1) = 1;    
    Bi = zeros(1,nvars,horzlim+1);
    if ii == target
    for l=0:horzlim
        C=M^l;
        if l>0
            Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars)+Bi(:,:,l); % cumulate because differences
        else
            Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars);
        end
    end
    else
    for l=0:horzlim
        C=M^l;
        Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars); % cumulate because differences
    end
    end
    %compute V initialize           
    FEVD(ii,:,1) = diag((Bi(:,:,1)*Atilde)'*(Bi(:,:,1)*Atilde))'; % 

    %impact of all shocks on error variance at each horizon
    for jj=1:horzlim-1 
        FEVD(ii,:,jj+1) = FEVD(ii,:,jj) + diag((Bi(:,:,jj+1)*Atilde)'*(Bi(:,:,jj+1)*Atilde))';

    end
    
    %Sum across rows should add to FEVD - divide each row at each horizon by
    %sum
    
end
FEVDcont = FEVD./sum(FEVD,2);
end
