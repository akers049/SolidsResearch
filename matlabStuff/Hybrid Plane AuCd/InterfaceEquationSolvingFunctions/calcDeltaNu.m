function [delta, nu] = calcDeltaNu(a, nHat, UI)
% This function just calculates the delta and nu given in euation 7.13 on
% page 113. It needs the a and nHat vector that satisfy the twinning
% equation, and the transformation matrix UI.

delta = dot(a, (UI*inv(UI^2 - eye(3))*nHat) );
nu = trace(UI^2) - det(UI^2) - 2 + sum(a.^2)/(2*delta);
end