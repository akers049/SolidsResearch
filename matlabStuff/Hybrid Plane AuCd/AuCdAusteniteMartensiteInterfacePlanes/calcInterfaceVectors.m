function [b, mHat] = calcInterfaceVectors(lamdas, eHats, kappa)
% This function calculates the austenite-martensite interface vectors b and
% mHat, given the eigenvalues (lamdas) and eigenvectors (eHats) of the C
% matrix that is calculated in equation 7.16 on page 113. It also needs the
% kappa value (-1 or 1)

eHat1 = eHats(:, 1);
eHat3 = eHats(:, 3);

ae1factor = sqrt(lamdas(3)*(1-lamdas(1))/(lamdas(3) - lamdas(1)));
ae3factor = sqrt((lamdas(1)*(lamdas(3)-1))/(lamdas(3)-lamdas(1)));

mHatfrontfactor = (sqrt(lamdas(3))-sqrt(lamdas(1)))/sqrt(lamdas(3)-lamdas(1));
mHate1factor = -sqrt(1-lamdas(1));
mHate3factor = sqrt(lamdas(3)-1);

m = mHatfrontfactor*(mHate1factor*eHat1 + kappa*mHate3factor*eHat3);

rho = norm(m);

mHat = m/rho;
  
b = rho*(ae1factor*eHat1 + kappa*ae3factor*eHat3);
end