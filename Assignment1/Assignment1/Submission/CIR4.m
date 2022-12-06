% Active Filte Circuit
% Output at node 4

% global   G C b


InitiateCircuit
r1 = 1200;
r2 = 1010;
r3 = 1140;
r4 = 1450;
r5 = 1540;

% vol('V1','1','0','type','AC','AMPLITUDE',1);
vol('V1','1','0','type','DC','VAL_DC',1);
res('r1','1','2',r1);
cap('C1','2','0',200e-9);
res('r2','2','0',r2);
vccs('G1','3','0','2','0',1/r3);
res('r4','3','0',r4);
vcvs('E1','4','0','3','0',-10);
res('r5','4','0',r5);


% InitiateCircuit
% vol('V1','n1','0','type','AC','AMPLITUDE',1)
% res('r1','n1','n2',100);
% res('r2','n2','0',50);