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

i = 1; 
delta_x = zeros(size(Xguess));%get matrix of deltax with random number as same size of Xguess

while 1
                        
    f = makeFvect(Xguess);
    J = make_nlJacobian(Xguess);
    phi_x = G * Xguess + f - alpha * Bdc;   %calculate phi using equations from slides
    Dphi_x = G + J;               %calculate d phi using G and J (jacobian matrix)
    delta_x = -Dphi_x\phi_x;      %calculate deltax using the equation from slides 
    dX(i) = norm(delta_x);
    Xguess = Xguess + delta_x;    %new guess
    Xdc = Xguess;  
    i = i + 1; 
    if norm(delta_x) < maxerr
        break
    end
end