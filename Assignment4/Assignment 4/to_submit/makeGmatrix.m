function Gmat = makeGmatrix()

global elementList nodeMap
Gmat = sparse(elementList.n,elementList.n);

%% add stamp for resistors

for I=1:elementList.Resistors.numElements
    
    % access nodes Numbers of the Resistor
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.Resistors.nodeNumbers(I,:);
    
    % get the conductance for the resisitor
    % resistance is stored in the field named value
    g = 1/( elementList.Resistors.value(I));
    
    if(nodes(1)~=0) && (nodes(2)~=0)
        Gmat(nodes(1),nodes(1)) = Gmat(nodes(1),nodes(1))  + g;
        Gmat(nodes(1),nodes(2)) = Gmat(nodes(1),nodes(2)) - g;
        Gmat(nodes(2),nodes(1)) = Gmat(nodes(2),nodes(1)) - g;
        Gmat(nodes(2),nodes(2)) = Gmat(nodes(2),nodes(2)) + g;
    elseif (nodes(1)==0) && (nodes(2)~=0)
        
        Gmat(nodes(2),nodes(2)) = Gmat(nodes(2),nodes(2)) + g;
        
    elseif (nodes(1)~=0) && (nodes(2)==0)
        Gmat(nodes(1),nodes(1)) = Gmat(nodes(1),nodes(1))  + g;
    end
    
end

%% for  DC voltage sources
for I=1:elementList.DC_VolSources.numElements
    % access nodes Numbers of the Resistor
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.DC_VolSources.nodeNumbers(I,:);
    
    % get the extra current indices
    nX = elementList.DC_VolSources.curIdx(I);
    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1))  + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) +1;
    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2))  - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end
    
end

%% for  AC voltage sources
for I=1:elementList.AC_VolSources.numElements
    % access nodes Numbers of the Resistor
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.AC_VolSources.nodeNumbers(I,:);
    
    % get the extra current indices
    nX = elementList.AC_VolSources.curIdx(I);
    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1))  + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) +1;
    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2))  - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end
    
end


%% for  TR voltage sources
for I=1:elementList.TR_VolSources.numElements
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.TR_VolSources.nodeNumbers(I,:);
    
    % get the extra current indices
    nX = elementList.TR_VolSources.curIdx(I);
    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1))  + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) +1;
    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2))  - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end
    
end
%% for  HB voltage sources
for I=1:elementList.HB_VolSources.numElements
    %nodes of a I^{th} element are located in Row I of the nodeNumbers
    %field
    
    nodes = elementList.HB_VolSources.nodeNumbers(I,:);
    
    % get the extra current indices
    nX = elementList.HB_VolSources.curIdx(I);
    if(nodes(1)~=0)
        Gmat(nX,nodes(1)) = Gmat(nX,nodes(1))  + 1;
        Gmat(nodes(1),nX) = Gmat(nodes(1),nX) +1;
    end
    if(nodes(2)~=0)
        Gmat(nX,nodes(2)) = Gmat(nX,nodes(2))  - 1;
        Gmat(nodes(2),nX) = Gmat(nodes(2),nX) - 1;
    end
    
end
%% Add for  Inductors  voltage sources
for I=1:elementList.Inductors.numElements
    nX= elementList.Inductors.curIdx(I);
    % get nodes
    n1 = elementList.Inductors.nodeNumbers(I,1);
    n2 = elementList.Inductors.nodeNumbers(I,2);
    
    if(n1 ~=0)
        Gmat(nX,nP) =Gmat(nX,nP)+1;
        Gmat(nP,nX) =Gmat(nP,nX)+1;
    end
    if(n2 ~=0)
        Gmat(nX,nN) =Gmat(nX,nN)-1;
        Gmat(nN,nX) =Gmat(nN,nX)-1;
    end
end

%% Add for Opamps
for I=1:elementList.OpAmps.numElements
    % get the current Index
    nX= elementList.OpAmps.curIdx(I);
    % get nodes
    nP = elementList.OpAmps.nodeNumbers(I,2); % positive node
    nN = elementList.OpAmps.nodeNumbers(I,3); % negative node
    nO = elementList.OpAmps.nodeNumbers(I,1); % Output node
    
% %         if (nN~=0) & (n3~=0)
% %             Gmat(nX,nN) = -1;
% %             Gmat(nX, n3) = 1;
% %         end
% %         if (n2~=0) & (n3==0)
% %             Gmat(nX, nN) = -1;
% %         end
% %         if (n2==0) & (n3~=0)
% %             Gmat(nX, n3) = -1;
% %         end
% %         if (n1~=0)
% %             Gmat(n1, nX) = 1;
% %         end

%     inplus =nP;inminus =nN; out=nO; xr = nX; Gain = 500;
%     if(inplus ~=0 && inplus ~=0)
%         Gmat(xr,inplus) = Gmat(xr,inplus)- Gain;
%         Gmat(xr,inminus) = Gmat(xr,inminus)+ Gain;
%         Gmat(xr,out) = Gmat(xr,out)+1;
%         Gmat(out,xr) =Gmat(out,xr) +1; %column
%     end
%     if(inplus~=0 && inminus==0)
%         Gmat(xr,inplus) = Gmat(xr,inplus)- Gain;
%         Gmat(xr,out) = Gmat(xr,out)+1;
%         Gmat(out,xr) =Gmat(out,xr) +1;%column
%     end
%     if(inplus==0 && inminus ~=0)
%         Gmat(xr,inminus) = Gmat(xr,inminus)+ Gain;
%         Gmat(xr,out) = Gmat(xr,out)+1;
%         Gmat(out,xr) =Gmat(out,xr) +1;%column
%     end
        if(nP ~=0)
            Gmat(nX,nP) = Gmat(nX,nP) -1;
        end
        if(nN ~=0)
            Gmat(nX,nN) = Gmat(nX,nN) +1;
        end
        if(nO ~=0)
            Gmat(nO,nX) =Gmat(nO,nX) +1;
        end
end
end