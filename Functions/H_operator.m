function mag = H_operator(psi1,psi2,dA,eps)


mag=psi1(:,:,1)'.*(dA.*eps.*psi2(:,:,1))+psi1(:,:,2)'.*(da.*psi2(:,:,2));


end
