%% NETLIST
InitiateCircuit;
vol('v1','1','0','Type', 'HB','AMPLITUDE',1, 'frequency',60)
res('R1','1','2',1e3)
diode('D1','2','3','Is',1e-14,'Vt',0.025)
res('R2','3','0',5000)
cap('C1','3','0',1e-6);

add_extraIndices;