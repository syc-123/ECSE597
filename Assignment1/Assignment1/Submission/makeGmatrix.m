function Gmat = makeGmatrix()

global elementList nodeMap
elementList.n = elementList.numNodes; % MNA size
Gmat = sparse(elementList.n,elementList.n);

% add stamp for resistors

for I=1:elementList.Resistors.numElements
    
    % access nodes Numbers of the Resistor
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.Resistors.nodeNumbers(I,:);
    
    
    % get the conductance for the resisitor
    % resistance is stored in the field named value
    g = 1/(elementList.Resistors.value(I));
    
    if(nodes(1)~=0) && (nodes(2)~=0)
        Gmat(nodes(1),nodes(1)) = Gmat(nodes(1),nodes(1)) + g;
        Gmat(nodes(1),nodes(2)) = Gmat(nodes(1),nodes(2)) - g;
        Gmat(nodes(2),nodes(1)) = Gmat(nodes(2),nodes(1)) - g;
        Gmat(nodes(2),nodes(2)) = Gmat(nodes(2),nodes(2)) + g;
    elseif (nodes(1)==0) && (nodes(2)~=0)
        
        Gmat(nodes(2),nodes(2)) = Gmat(nodes(2),nodes(2)) + g;
       
    elseif (nodes(1)~=0) && (nodes(2)==0)
        Gmat(nodes(1),nodes(1)) = Gmat(nodes(1),nodes(1)) + g;
    end
    
end

% for  DC voltage sources

for I=1:elementList.DC_VolSources.numElements
    % access nodes Numbers of the Resistor
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.DC_VolSources.nodeNumbers(I,:);
    
    % For voltage sources we need to add an extra row and a column 
    nX= elementList.n+1;  % extra node 
    elementList.n = elementList.n+1;  % update n
        
    % store the extra node as you will need this in in filling Bvector
    % to do this we use the field currIndex in DC_volSources
    elementList.DC_VolSources.currIndex(I) = nX;
    Gmat(nX,nX) = 0;

    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1)) + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) + 1;

    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2)) - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end 

    
end 


% add stamps for other elements here...

% for inductors 
for I=1:elementList.Inductors.numElements
    nodes = elementList.Inductors.nodeNumbers(I,:);
    nX= elementList.n+1;
    elementList.n= elementList.n+1;
    elementList.Inductors.currIndex(I) = nX;
    Gmat(nX,nX)=0;
    
    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1)) + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) + 1;

    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2)) - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end 
end

% for voltage controlled voltage source
for I=1:elementList.VCVS.numElements
    nodes = elementList.VCVS.nodeNumbers(I,:);
    nX= elementList.n+1;
    elementList.n= elementList.n+1;
    A = elementList.VCVS.value(I);
    Gmat(nX,nX)=0;
    
    if(nodes(3)~=0)
        Gmat(nX,nodes(3)) = Gmat(nX,nodes(3)) - A;
    end
    if(nodes(4)~=0)
        Gmat(nX,nodes(4)) = Gmat(nX,nodes(4)) + A;
    end
    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1)) + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) + 1;
    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2)) - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end
end

% for voltage controlled current source
for I=1:elementList.VCCS.numElements
    nodes = elementList.VCCS.nodeNumbers(I,:);
    gm = elementList.VCCS.value(I);

    if(nodes(1)~=0)
        if(nodes(3)~=0)
            Gmat(nodes(1),nodes(3)) = Gmat(nodes(1),nodes(3)) + gm;
        end
        if(nodes(4)~=0)
            Gmat(nodes(1),nodes(4)) = Gmat(nodes(1),nodes(4)) - gm;
        end
    end

    if(nodes(2)~=0)
        if(nodes(3)~=0)
            Gmat(nodes(2),nodes(3)) = Gmat(nodes(2),nodes(3)) - gm;
        end
        if(nodes(4)~=0)
            Gmat(nodes(2),nodes(4)) = Gmat(nodes(2),nodes(4)) + gm;
        end
    end
end

% for current controlled voltage source
for I=1:elementList.CCVS.numElements
    nodes = elementList.CCVS.nodeNumbers(I,:);
    nX= elementList.n+1;
    elementList.n= elementList.n+1;
    L = elementList.CCVS.value(I);
    Gmat(nX,nX)=0;

    Gmat(nX,nodes(1)) = Gmat(nX,nodes(1)) + 1;
    Gmat(nX,nodes(2)) = Gmat(nX,nodes(2)) - 1;
    Gmat(nodes(1),nx) = Gmat(nodes(1),nx) + 1;
    Gmat(nodes(2),nx) = Gmat(nodes(2),nx) - 1;

    nX= elementList.n+1;
    elementList.n= elementList.n+1;
    Gmat(nX,nX)=0;

    Gmat(nx,nodes(3)) = Gmat(nx,nodes(3)) + 1; 
    Gmat(nx,nodes(4)) = Gmat(nx,nodes(4)) - 1; 
    Gmat(nodes(3),nx) = Gmat(nodes(3),nx) + 1; 
    Gmat(nodes(4),nx) = Gmat(nodes(4),nx) - 1;
    Gmat(nx,nx-1) = Gmat(nx,nx-1)-L
    
end

% for current controlled current source
for I=1:elementList.CCCS.numElements
    nodes = elementList.CCCS.nodeNumbers(I,:);
    nX= elementList.n+1;
    elementList.n= elementList.n+1;
    A = elementList.CCCS.value(I);
    Gmat(nX,nX)=0;

    Gmat(nX,nodes(1)) = Gmat(nX,nodes(1)) + 1;
    Gmat(nX,nodes(2)) = Gmat(nX,nodes(2)) - 1;
    Gmat(nodes(1),nx) = Gmat(nodes(1),nx) + 1; 
    Gmat(nodes(2),nx) = Gmat(nodes(2),nx) - 1; 
    Gmat(nodes(3),nx) = Gmat(nodes(3),nx) + A; 
    Gmat(nodes(4),nx) = Gmat(nodes(4),nx) - A; 

end

% % for mutual inductor
% for I=1:elementList.Mutual_Ind.numElements
%     nodes = elementList.Mutual_Ind.nodeNumbers(I,:);
%     nX= elementList.n+1;
%     elementList.n= elementList.n+1;
%     Gmat(nX,nX)=0;
% 
%     Gmat(nX,nodes(1)) = Gmat(nX,nodes(1)) + 1;
%     Gmat(nX,nodes(2)) = Gmat(nX,nodes(2)) - 1;
%     Gmat(nodes(1),nx) = Gmat(nodes(1),nx) + 1;
%     Gmat(nodes(2),nx) = Gmat(nodes(2),nx) - 1;
%     
%     nX= elementList.n+1;
%     elementList.n= elementList.n+1;
%     Gmat(nX,nX)=0;
%     elementList.Mutual_Ind.currIndex(I)=nX
% 
%     Gmat(nx,nodes(3)) = Gmat(nx,nodes(3)) + 1; 
%     Gmat(nx,nodes(4)) = Gmat(nx,nodes(4)) - 1; 
%     Gmat(nodes(3),nx) = Gmat(nodes(3),nx) + 1; 
%     Gmat(nodes(4),nx) = Gmat(nodes(4),nx) - 1;
%         
% end
% 
% for ideal opamp
for I=1:elementList.OpAmps.numElements
    nodes = elementList.OpAmps.nodeNumbers(I,:);
    nX = elementList.n+1;
    elementList.n = elementList.n+1;
    Gmat(nX,nX)=0;
    
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2))-1;
    end
    if(nodes(3)~=0)
        Gmat(nX,nodes(3)) = Gmat(nX,nodes(3))+1;
    end
    if(nodes(1)~=0)
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX)+1;
    end

end


end 