function Cmat = makeCmatrix()

global elementList nodeMap
n = elementList.n
Cmat = sparse(n,n);

% add stamps here...

% for Capacitors  
for I=1:elementList.Capacitors.numElements
    nodes = elementList.Capacitors.nodeNumbers(I,:);
    c = elementList.Capacitors.value(I);

    if(nodes(1)~=0) && (nodes(2)~=0)
        Cmat(nodes(1),nodes(1)) = Cmat(nodes(1),nodes(1)) + c;
        Cmat(nodes(1),nodes(2)) = Cmat(nodes(1),nodes(2)) - c;
        Cmat(nodes(2),nodes(1)) = Cmat(nodes(2),nodes(1)) - c;
        Cmat(nodes(2),nodes(2)) = Cmat(nodes(2),nodes(2)) + c;
    elseif (nodes(1)==0) && (nodes(2)~=0)        
        Cmat(nodes(2),nodes(2)) = Cmat(nodes(2),nodes(2)) + c;       
    elseif (nodes(1)~=0) && (nodes(2)==0)
        Cmat(nodes(1),nodes(1)) = Cmat(nodes(1),nodes(1)) + c;
    end

end

% for inductors 
for I=1:elementList.Inductors.numElements
    idx = elementList.Inductors.currIndex(I)
    l = elementList.Inductors.value(I);
    Cmat(idx,idx)= -l;
    
end

% % for mutual inductor
% for I=1:elementList.mutual_ind.numElements
%     idx = elementList.mutual_ind.currIndex(I)
%     Cmat(idx,idx) = 0;
%     L = elementList.mutual_ind.value(I)
%     nodes = elementList.mutual_ind.nodeNumbers(I,:);
% 
%     Cmat(idx-1,idx-1) = Cmat(idx-1,idx-1) - L(1);
%     Cmat(idx-1,idx) = Cmat(idx-1,idx) - L(3);
%     Cmat(idx,idx-1) = Cmat(idx,idx-1) - L(3);
%     Cmat(idx,idx) = Cmat(idx,idx) - L(2);
%       
% end

% add stamp for other elelemnts here 