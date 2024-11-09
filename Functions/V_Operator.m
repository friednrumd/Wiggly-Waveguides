function mag = V_Operator(psi1,psi2,dA,eps)

mag= 1/2*trace(       psi1(:,1:2,1)'  *(dA.*eps.*psi2(:,1:2,1))     + psi1(:,1:2,2)'*(dA.*psi2(:,1:2,2)) );

end