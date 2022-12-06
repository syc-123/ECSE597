function Bvec = makeBvector()

global elementList 

numNodes = elementList.n;
Bvec = sparse(numNodes,1);


% add stamps here...
for I=1:elementList.DC_VolSources.numElements
    % acess the current  index 
    idx = elementList.DC_VolSources.currIndex(I);
    Bvec(idx, 1) = elementList.DC_VolSources.VAL_DC(I);
end

for I=1:elementList.DC_CurSources.numElements
    nodes = elementList.DC_CurSources.nodeNumbers(I,:);
    cur = elementList.DC_CurSources.value(I);

    if(nodes(1)~=0)
        Bvec(nodes(1),1) = Bvec(nodes(1),1) + cur;
    end
    if(nodes(2)~=0)
        Bvec(nodes(2),1) = Bvec(nodes(2),1) - cur;
    end 
end

end 
