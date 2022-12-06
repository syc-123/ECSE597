%% NETLIST
InitiateCircuit

vol('V1','1','0','type','DC', 'VAL_DC',10)
res('R1','1','2',4.7e3)
res('R2','4','0',3.3e3)
% BJT('Name',CollectorNode,baseNode,EmitterNode, 'type', 'NPN', 'IS', 1e-15,'Vt', 0.025,'alphaF',0.99,'alphaR',0.1)
bjt('BJT1','2','3','4', 'type', 'NPN', 'IS', 1e-15,'Vt', 0.025,'alphaF',0.99,'alphaR',0.1)


%----------------------------------------------------------------------
%----------------------------------------------------------------------
% Base to Ground %4V Active % 6V Saturation %0 Cut-off
Vb = 4; 

vol('V2','3','0','type','DC', 'VAL_DC',Vb)

%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% create MNA.
% add_extraIndices
G= makeGmatrix;

%% Compare the performance of dcsolve vs dcsolvecont
Xguess = zeros(size(G,1),1);
maxerr = 1e-4;

[Xdcsolve, dX] = dcsolve(Xguess,maxerr)
 
 
 Xdcsolve_cont = dcsolvecont(50,maxerr)
