function r = fsolve( fpoints ,out)
%  fsolve(fpoints ,out)
%  Obtain frequency domain response

% fpoints: are the frequency points
%out is the output node name 
out_NodeNumber = getNodeNumber(out) ;  % node number associated with each node 

Gmat = makeGmatrix;
Cmat = makeCmatrix;
[~, Bac] = makeBvector;
for I=1:length(fpoints)
    
    A = Gmat + 2*pi*fpoints(I)*1i*Cmat;
    
    x = A\Bac;
    r(I,1) = x(out_NodeNumber);
end 
    
    
      