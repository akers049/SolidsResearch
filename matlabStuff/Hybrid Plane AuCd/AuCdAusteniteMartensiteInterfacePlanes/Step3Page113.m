function [bs, mHats] = Step3Page113(lam, a, nHat, UI)
% This function executes step 3 on page 113, given the lamda value
% (lam), the a and nHat vectors that satisfy the twinning equation, and the
% transformation matrix UI.
% It outputs cells bs and mHats that contain the solution vectors b and
% mHat, respectively to the austenite martensite twinning equation.

C1 = (UI + lam*(nHat*transpose(a)))*(UI + lam*(a*transpose(nHat)));
[lamdas, eHats] = calcAndOrderEigs(C1);
[b1, mHat1] = calcInterfaceVectors(lamdas, eHats, 1);
[b2, mHat2] = calcInterfaceVectors(lamdas, eHats, -1);

bs = {b1; b2};
mHats = {mHat1; mHat2};
end