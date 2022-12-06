function Fvect = makeFvect(X)
% this function creates a nonlinear vector in the MNA matrices.
% X is the MNa vector containing the node voltages/currents/flux/charge

global elementList

Fvect  = sparse(elementList.n,1);

% add diodes using this
for I=1:elementList.Diodes.numElements
    nodes = elementList.Diodes.nodeNumbers(I,:);
    % nodes(1) is the positive node of the diode
    % node(2) is the negative node of the diode.
    
    n1 = nodes(I,1);
    n2 = nodes(I,2);
   
     % get other parameters
     Vt = elementList.Diodes.Vt(I);  % thermal voltage of the I^th Diode 
     Is = elementList.Diodes.Is(I);       % sturation current for the I^th Diode
     
    if(nodes(1)~=0 && nodes(2)~=0)
        Vd = X(n1)-X(n2);
        
        Fvect(n1,1) = Fvect(n1,1) + Is*( exp(Vd/Vt) - 1) ;
        Fvect(n2,1) = Fvect(n2,1) - Is*( exp( Vd/Vt) - 1) ;
    elseif (nodes(1)==0 && nodes(2)~=0)
        
        Vd = -X(n2);
        Fvect(n2,1) = Fvect(n2,1) - Is*( exp( Vd/Vt) - 1) ;
        
    elseif(nodes(1)~=0 && nodes(2)==0)
        Vd = X(n1);  
        Fvect(n1,1) = Fvect(n1,1) + Is*( exp(Vd/Vt) - 1) ;
    end
    
end

% add BJts using this
for I=1:elementList.BJTs.numElements
   
     
     %get Nodes Numbers
     cNode = elementList.BJTs.nodeNumbers(I,1);  % collector node.
     bNode = elementList.BJTs.nodeNumbers(I,2);    % base node.
     eNode = elementList.BJTs.nodeNumbers(I,3);    % emitter node
     
     % get other parameters
     Vt = elementList.BJTs.Vt(I);  % thermal voltage of the I^th BJT
     Is = elementList.BJTs.Is(I);       % sturation current for the I^th BJT
     alphaR = elementList.BJTs.alphaR(I);
     alphaF =  elementList.BJTs.alphaF(I);
     
 
    
     if(cNode~=0)&&(bNode~=0)&&(eNode~=0) % all nodes present
         % get nodal voltages
         Vbe = X(bNode) -X(eNode);
         Vbc = X(bNode) -X(cNode);
         % diode currents 
         If = Is*(exp(Vbe/Vt)-1);
         Ir = Is*(exp(Vbc/Vt)-1);
         
      Fvect(cNode) =  Fvect(cNode)  -Ir            +    alphaF*If ;
      Fvect(bNode) =  Fvect(bNode)  +Ir*(1-alphaR) +    If*(1-alphaF);
      Fvect(eNode) =  Fvect(eNode)  +Ir*alphaR     -    If;
     elseif (cNode~=0)&&(bNode==0)&&(eNode~=0) % Base is Grounded
 % get nodal voltages
         Vbe =  -X(eNode);
         Vbc =  -X(cNode);
         % diode currents 
         If = Is*(exp(Vbe/Vt)-1);
         Ir = Is*(exp(Vbc/Vt)-1);
         
      Fvect(cNode) =  Fvect(cNode)  -Ir            +    alphaF*If ;
      Fvect(eNode) =  Fvect(eNode)  +Ir*alphaR     -    If;
    
     elseif (cNode~=0)&&(bNode~=0)&&(eNode==0) % Emitter is Grounded
           % get nodal voltages
         Vbe = X(bNode) ;
         Vbc = X(bNode) -X(cNode);
         % diode currents 
         If = Is*(exp(Vbe/Vt)-1);
         Ir = Is*(exp(Vbc/Vt)-1);
         
      Fvect(cNode) =  Fvect(cNode)  -Ir            +    alphaF*If ;
      Fvect(bNode) =  Fvect(bNode)  +Ir*(1-alphaR) +    If*(1-alphaF);
     
      
    elseif (cNode==0)&&(bNode~=0)&&(eNode~=0) % Collector is Grounded
       % get nodal voltages
         Vbe = X(bNode) -X(eNode);
         Vbc = X(bNode);
         % diode currents 
         If = Is*(exp(Vbe/Vt)-1);
         Ir = Is*(exp(Vbc/Vt)-1);
         
     
      Fvect(bNode) =  Fvect(bNode)  +Ir*(1-alphaR) +    If*(1-alphaF);
      Fvect(eNode) =  Fvect(eNode)  +Ir*alphaR     -    If;
      
     elseif(cNode~=0)&&(bNode==0)&&(eNode==0) % Base and  Emitter are grounded
        % fill this up  
     end 
end