function mag = conjS_Metric(psi1,psi2,dA)
    
mag=1/2*( conj(psi1(:,1,1))'*( dA.*psi2(:,2,2) )-conj(psi1(:,2,1))'*( dA.*psi2(:,1,2) ) )-1/2*( conj(psi1(:,1,2))'*( dA.*psi2(:,2,1) )-conj(psi1(:,2,2))'*( dA.*psi2(:,1,1) ) );


end



