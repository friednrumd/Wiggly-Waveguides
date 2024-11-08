function mag = H_Operator(psi1,psi2,dA,eps)

mag=sum( psi1(:,:,1)'*(dA.*eps.*psi2(:,:,1))+psi1(:,:,2)'*(dA.*eps.*psi2(:,:,2))   ,"All");

end