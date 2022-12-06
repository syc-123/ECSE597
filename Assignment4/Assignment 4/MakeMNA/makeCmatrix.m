function Cmat = makeCmatrix()

global elementList nodeMap

n = elementList.n;
Cmat = sparse(n,n);




% add stamp for resistors

for I=1:elementList.Capacitors.numElements
    
    % access nodes Numbers of the Resistor
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.Capacitors.nodeNumbers(I,:);

    % get the conductance for the resisitor
    % resistance is stored in the field named value
   c= elementList.Capacitors.value(I);
    
    if(nodes(1)~=0) && (nodes(2)~=0)
        Cmat(nodes(1),nodes(1)) = Cmat(nodes(1),nodes(1))  + c;
        Cmat(nodes(1),nodes(2)) = Cmat(nodes(1),nodes(2)) - c;
        Cmat(nodes(2),nodes(1)) = Cmat(nodes(2),nodes(1)) - c;
        Cmat(nodes(2),nodes(2)) = Cmat(nodes(2),nodes(2)) + c;
    elseif (nodes(1)==0) && (nodes(2)~=0)
        
        Cmat(nodes(2),nodes(2)) = Cmat(nodes(2),nodes(2)) + c;
        
    elseif (nodes(1)~=0) && (nodes(2)==0)
        Cmat(nodes(1),nodes(1)) = Cmat(nodes(1),nodes(1))  + c;
    end
    
end


% for inductors 
for I=1:elementList.Inductors.numElements
    
end

% add stamp for other elelemnts here 