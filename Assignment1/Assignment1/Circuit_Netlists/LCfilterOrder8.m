% LC filter 8th order
% output at node 6

 InitiateCircuit

R1 = 50;	
L1 = 0.00310511;	
C1 = 3.53642e-06;
L2 = 0.0132337;	
C2 = 6.24524e-06;	
L3 = 0.0156131;
C3 = 5.29349e-06; 	
L4 = 0.00884106;	
C4 = 1.24205e-06;	
R2 = 50;

vol('V1','1','0','type','AC','AMPLITUDE',1)

res('R1','1','2',R1)
ind('L1','2','3',L1)
cap('C1','3','0',C1)
ind('L2','3','4',L2)
cap('C2','4','0',C2)
ind('L3','4','5',L3)
cap('C3', '5','0',C3)
ind('L4','5','6',L4)
cap('C4','6','0',C4)
res('R2','6','0',R2)