function r = nonlinear_fsolve(fpoints ,out)

% This function  Obtains frequency response of the nonlinear function
% global variables G C b bac
% Inputs: 1. fpoints is a vector containing the fequency points\
%         2. out is the node name of the output node. 
             
% Outputs: r is a vector containing the value of
%            of the response at the points fpoint

global  elementList

outNodeNumber = getNodeNumber(out);


%compute the Jacoba
G = makeGmatrix;
C = makeCmatrix;
[Bdc, Bac] = makeBmatrix;
J =  G+ make_nlJacobian(Xdc);


for I = 1: length(fpoints)
    Xr = (J+2*pi*1i*fpoints(I)*C)\Bac;
    r(I) = Xr(outNodeNumber);
end