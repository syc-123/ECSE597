InitiateCircuit

vol('V1', 'n1', 'gnd', 'type', 'AC', 'AMPLITUDE', 2)
cap('C1', 'n1', 'n2', 220e-9);
cap('C2', 'n2', 'n3', 220e-9);
res('R2', 'n3', 'gnd', 10e3);
res('R1', 'n2', 'n4', 10e3);
res('Rout', 'n4', 'gnd',500);
opamp('E1', 'n3', 'n4','n4');