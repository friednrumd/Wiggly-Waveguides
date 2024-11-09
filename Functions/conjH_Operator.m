function mag = conjH_Operator(psi1,psi2,dA,eps)

mag= 1/2*trace( conj( psi1(:,:,1) )'*(dA.*eps.*psi2(:,:,1))+conj( psi1(:,:,2) )'*( dA.*psi2(:,:,2)) );


end
