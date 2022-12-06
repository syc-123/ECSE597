%% NeTLIST
InitiateCircuit

vol('v1','1','0','type','dc','val_dc',1)

diode('D1','1','2','IS',1e-15,'Vt',26e-3);
res('r1','2','0',50)

% add_extraIndices
G = makeGmatrix;
n = elementList.n;% MNA size 




%%  1.   Call your  DC solve function  here and report your observations in your report
Xguess = zeros(n,1);
maxerr = 1e-4;
[Xdc, dX] = dcsolve(Xguess,maxerr)