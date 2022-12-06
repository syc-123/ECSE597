function [tpoints, y]= FE_method(tEnd,h, outNode)
% This function uses FORWARD EULER method to compute the transient reponse
% of the circuit.
%Inputs: 1. tEnd:  The simulation starts at time = 0s  and ends at 
%                  time = tEND s.
%        2. h: length  of step size.
%        3. outNode: is the node for which the transient is required.
%Output:  1.y:  is the transient response at outNode.
%
%
%Note: The function stub provided above is just an example. You can modify the 
%      in function in any fashion. 
%--------------------------------------------------------------------------

global elementList

out_NodeNumber = getNodeNumber(outNode) ;
tpoints = 0:h:tEnd; 

Gmat = makeGmatrix;
Cmat = makeCmatrix;

x_n = zeros(size(Gmat,2),1);
b_n = zeros(size(Gmat,1),1);
y = zeros(1,length(tpoints));


for I=1:length(tpoints)
    b_n1 = makeBt(tpoints(I));
    x_n1 = (Cmat) \ (h * (b_n - Gmat * x_n)) + x_n;
    y(I) = x_n1(out_NodeNumber);
    x_n = x_n1;
    b_n = b_n1;
    
end 
     
    
end     


