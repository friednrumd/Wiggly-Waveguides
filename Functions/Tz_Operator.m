function mag = Tz_Operator(psi1,psi2,dA,eps)

mag= 1/2*(       psi1(:,3,1)'  *(dA.*eps.*psi2(:,3,1))     + psi1(:,3,2)'*(dA.*psi2(:,3,2)) );

end