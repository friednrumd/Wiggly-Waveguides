function df = conjS_Metric(psi1,psi2,dA)
    
mag=0;


mag=mag+1/2*conj(psi1(:,1,1))'*( dA.*psi2(:,2,2) );
mag=mag-1/2*conj(psi1(:,2,1))'*( dA.*psi2(:,1,2) );


mag=mag-1/2*conj(psi1(:,1,2))'*( dA.*psi2(:,2,1) );
mag=mag+1/2*conj(psi1(:,2,2))'*( dA.*psi2(:,1,1) );

df=mag;

end


