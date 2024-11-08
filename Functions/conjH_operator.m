function mag = conjH_operator(psi1,psi2,dA,eps)


mag=conj(psi1(:,:,1))'.*(dA.*eps.*psi2(:,:,1))+conj(psi1(:,:,2))'.*(da.*psi2(:,:,2));


end
