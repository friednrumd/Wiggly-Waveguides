function q = zIndIterator(Hamiltonian,vec0,eig0,r0,dw)

    [vec1, eig1]=Perturber1(Hamiltonian*dw/2,vec0,eig0);
    [vec2, eig2]=Perturber1(Hamiltonian*dw/2,vec0+vec1*dw/2,eig0+eig1*dw/2);
    [vec3, eig3]=Perturber1(Hamiltonian*dw,vec2,eig2);

    [vecf, eigf]=([vec3, eig3]+[vec2, eig2]+[vec1, eig1]+)

end

