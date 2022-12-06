function [Xdc dX] = dcsolve(Xguess,maxerr)
% Compute dc solution using newtwon iteration
% input: Xguess is the initial guess for the unknown vector. 
%        It should be the correct size of the unknown vector.
%        maxerr is the maximum allowed error. Set your code to exit the
%        newton iteration once the norm of DeltaX is less than maxerr
% Output: Xdc is the correction solution
%         dX is a vector containing the 2 norm of DeltaX used in the 
%         newton Iteration. the size of dX should be the same as the number
%         of Newton-Raphson iterations. See the help on the function 'norm'
%         in matlab. 
global elementList

[Bdc, Bac] = makeBvector;
G = makeGmatrix;

%maxerr = 1e-6;
I = 1;

Phi = G*Xguess + makeFvect(Xguess) - Bdc;
J = G + make_nlJacobian(Xguess);
deltaX = -J\Phi;
Xguess = Xguess + deltaX;
dX(I)=norm(deltaX);

while norm(deltaX)>maxerr
    f=makeFvect(Xguess);
    Phi = G*Xguess + f - Bdc;
    J = G + make_nlJacobian(Xguess);
    deltaX = -J\Phi;
    Xguess = Xguess + deltaX;
    I=I+1;
    dX(I) = norm(deltaX);
end

Xdc = Xguess;
