function r = fsolve(Gmat, Cmat, fpoints ,out)
%  fsolve(fpoints ,out)
%  Obtain frequency domain response

% fpoints: are the frequency points
%out is the output node name 
% out_NodeNumber = getNodeNumber({out}) ;  % node number associated with each node 

L = length(fpoints);
b = makeBmatrix;
f_response = zeros(1,L);

for i = 1:L
    
    A = Gmat + 2 * pi * fpoints(i) * 1i * Cmat;
    
    x = A \ b;
    
    f_response(1,i) = x(out,1);
    
end

r = f_response;