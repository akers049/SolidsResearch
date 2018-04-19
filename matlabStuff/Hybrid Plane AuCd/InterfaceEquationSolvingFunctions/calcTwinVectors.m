function [a, nHat] = calcTwinVectors(lamdas, eHats, G, kappa)
% This functuon calculates the vectors a and nHat that are the solutions to
% the twinning equation, given the eigenvalues (lamdas) and eigenvectors 
% (eHats) found in step 3 on page 69, that satisfy step 4. It also needs 
% the J transoformation matrix (G), and the kappa value (-1 or 1)

eHat1 = eHats(:, 1);
eHat3 = eHats(:, 3);
GT = transpose(G);

ae1factor = sqrt(lamdas(3)*(1-lamdas(1))/(lamdas(3) - lamdas(1)));
ae3factor = sqrt((lamdas(1)*(lamdas(3)-1))/(lamdas(3)-lamdas(1)));

nHatfrontfactor = (sqrt(lamdas(3))-sqrt(lamdas(1)))/sqrt(lamdas(3)-lamdas(1));
nHate1factor = -sqrt(1-lamdas(1))*GT;
nHate3factor = sqrt(lamdas(3)-1)*GT;
    
n = nHatfrontfactor*(nHate1factor*eHat1 + kappa*nHate3factor*eHat3);
rho = norm(n);
nHat = n/rho;

a = rho*(ae1factor*eHat1 + kappa*ae3factor*eHat3);
end