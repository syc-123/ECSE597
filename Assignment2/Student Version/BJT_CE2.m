InitiateCircuit

vol('v1','1','0','type','ac','amplitude',0.010)
res('Rs','1','2',50)% source resistor
cap('Cb','2','3',1e-6)% coupling capacitor 
res('r2','3','4',20e3)
vol('vcc','4','0','type','dc','val_dc',12)
res('r3','3','0',3.6e3 )
bjt('q1','5','3','6', 'type','NPN','Is',1e-13,'Vt',0.025,'alphaF',0.99,'alphaR',0.1)
res('r4','6','0',220)%emitter res
cap('Ce','6','0',100e-6)
res('Rc','5','4',1.2e3)% collector


cap('Cc','5','7',1e-6)

res('RL','7','0',10e3)% load 


% adding the Cmu and Cpi ( these are the parasitic Capacitances- Not shown in the diagram)
cap('Cmu','5','3',1e-9)
cap('Cpi','3','6',2e-9)
 

% add_extraIndices
G = makeGmatrix;
n = elementList.n;% MNA size 




%%  1.   Call your  DC solve function  here and report your observations in your report
% Xguess = zeros(n,1);
% [Xdc, dX] = dcsolve(Xguess,1e-3)  

%%  2.   Call your  Continutation  DC solve function.
Xdc = dcsolvecont(100,1e-6)


%%  3.   Test yout nonlinear frequency reponse function
fpoints = logspace(0,6,5000);
outNode ='7';


 % call your nonlinear frequnecy reposnse function here
        % Note: name the output of the nonlinear_fsolve as r    
r = nonlinear_fsolve(fpoints,outNode);
          
          

%  plot the results obtained from nonlinear_foslve for the outNode


close all
semilogx(fpoints, abs(r),'b')
grid on 
xlabel('Frequency');
ylabel('Magnitude of voltage at out node');


