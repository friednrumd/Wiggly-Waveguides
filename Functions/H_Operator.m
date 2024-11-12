function mag = H_Operator(psi1,psi2,dA,eps,d)

mag= 1/2*(       psi1(:,d,1).'  *(dA.*eps.*psi2(:,d,1))     + psi1(:,d,2).'*(dA.*psi2(:,d,2)) );

end