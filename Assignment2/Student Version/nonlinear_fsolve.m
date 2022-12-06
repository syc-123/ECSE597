function r = nonlinear_fsolve(fpoints ,out)

% This function  Obtains frequency response of the nonlinear function
% global variables G C b bac
% Inputs: 1. fpoints is a vector containing the fequency points\
%         2. out is the node name of the output node. 
             
% Outputs: r is a vector containing the value of
%            of the response at the points fpoint

global  elementList

outNodeNumber = getNodeNumber(out);



G = makeGmatrix;
C = makeCmatrix;
[Bdc, Bac] = makeBvector;
J = make_nlJacobian(dcsolvecont(100,1e-6));

L = length(fpoints);
f_response = zeros(1,L);

for i =1:L

    A = G + J + 2 * pi * fpoints(1,i) * 1i * C;
    x = A \ Bac 
    f_response(1,i) = x(outNodeNumber,1);

end

r = f_response;
