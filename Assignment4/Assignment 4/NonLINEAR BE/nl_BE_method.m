function [tpoints, X]= nl_BE_method(tEnd,h)
% This function uses BACKWARD EULER method to compute the transient reponse
% of a NONLINEAR circuit.
%Inputs: 1. tEnd:  The simulation starts at time = 0s  and ends at time =
%                  tEND s.
%        2. h: length  of step size.
%        3. outNode: is the node for which the transient is required.
%Output:  1. tpoints: list of time points.
%         2. y:  is the transient response at output node.
%
%Note: The function stub provided above is just an example. You can modify the
%      in function in any fashion.
%--------------------------------------------------------------------------

global elementList

maxerr = 1e-4;

%out_NodeNumber = getNodeNumber(outNode);
tpoints = 0:h:tEnd;

Gmat = makeGmatrix;
Cmat = makeCmatrix;

A = Gmat + Cmat / h;

x_n = zeros(size(Gmat,2),1);
X = zeros(size(Gmat,2),length(tpoints));

for I=1:length(tpoints)
    b_n1 = makeBt(tpoints(I));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xguess = x_n;
    flag = 1;

    while flag
        F = makeFvect(Xguess);
        J = make_nlJacobian(Xguess);
        
        phi_x = A * Xguess + F - (b_n1 + Cmat * x_n / h);
        Dphi_x = A + J;
        delta_x = -Dphi_x\phi_x;
        Xguess = Xguess + delta_x;
        x_n1 = Xguess;
        if norm(delta_x) < maxerr
            flag = 0;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X(I,:) = x_n1;
    x_n = x_n1;
    
end



end
