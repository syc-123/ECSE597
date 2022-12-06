InitiateCircuit

vol('V1', 'n1', 'gnd', 'type', 'AC', 'AMPLITUDE', 1)
res('R1', 'n1', 'n2', 100);
res('R2', 'n2', 'n3', 100);
cap('C2', 'n2', 'n4', 10e-6);
cap('C1', 'n3', 'gnd',10e-6);

opamp('E1', 'n5', 'n3','n4'); %opamp( Name, outNode, Positive Node, Neg Node)
res('R4', 'n4', 'n5', 50e3);
res('R3', 'n5', 'gnd', 5e3);
res('Rout', 'n4', 'gnd',1e3);