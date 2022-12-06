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
i = 1; 
delta_x = zeros(size(Xguess));%get matrix of deltax with random number as same size of Xguess

while 1
    f = makeFvect(Xguess);
    J = make_nlJacobian(Xguess);
    phi_x = G * Xguess + f - Bdc;
    Dphi_x = G + J;
    delta_x = -Dphi_x\phi_x;
    dX(i) = norm(delta_x);
    Xguess = Xguess + delta_x;
    Xdc = Xguess;
    i = i+1;
    if norm(delta_x) < maxerr
        break
    end

end

