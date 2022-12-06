function Xdc = dcsolvecont(n_steps,maxerr)
% Compute dc solution using newtwon iteration and continuation method
% (power ramping approach)
% inputs:
% n_steps is the number of continuation steps between zero and one that are
% to be taken. For the purposes of this assigments the steps should be 
% linearly spaced (the matlab function "linspace" may be useful).
% maxerr is the stopping criterion for newton iteration (stop iteration
% when norm(deltaX)<maxerr

global elementList

% size of MNA matrix
n = elementList.n;

i = 1;
alpha = linspace(0,1,n_steps);
Xguess = zeros(n,1);
 
while i <= n_steps
    Xguess = dcsolvealpha(Xguess,alpha(i),maxerr);
    Xdc = Xguess;
    i = i + 1;
end