function Xdc = dcsolvealpha(Xguess,alpha,maxerr)
% Compute dc solution using newtwon iteration for the augmented system
% G*X + f(X) = alpha*b
% Inputs: 
% Xguess is the initial guess for Newton Iteration
% alpha is a paramter (see definition in augmented system above)
% maxerr defined the stopping criterion from newton iteration: Stop the
% iteration when norm(deltaX)<maxerr
% Oupputs:
% Xdc is a vector containing the solution of the augmented system

global elementList 
[Bdc, Bac] = makeBvector;
G = makeGmatrix;

I = 1;
Phi = G*Xguess +  - alpha*Bdc;
J = G + make_nlJacobian(Xguess);
deltaX = -J\Phi;
Xguess = Xguess + deltaX;
dX(I) = norm(deltaX);

while norm(deltaX)>maxerr
    f =makeFvect(Xguess);
    Phi = G*Xguess + f - alpha*Bdc;
    J = G + make_nlJacobian(Xguess);
    deltaX = -J\Phi;
    Xguess = Xguess + deltaX;
    I=I+1;
    dX(I) = norm(deltaX);
end

Xdc = Xguess;