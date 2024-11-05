function [dvec,deig] = Perturber1(Hamiltonian,vecs,eigs)
   
deig=+Hamiltonian;
dvec=-(     vecs'*Hamiltonian*vecs    )/(       eigs-eigs'+diag(Inf*(1:length(eigs)))       );

end

