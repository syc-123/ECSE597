function J = make_nlJacobian(X)
% Input: X is the MNA solution vector. It contains the 
% node voltages and currents etc.

global elementList

J  = sparse(elementList.n,elementList.n);

% add diodes using this 
for I=1:elementList.Diodes.numElements
    nodes = elementList.Diodes.nodeNumbers(I,:);
    % nodes(1) is the positive node of the diode
    % node(2) is the negative node of the diode.
    n1 = nodes(1);
    n2 = nodes(2);

    % get other parameters
    Vt = elementList.Diodes.Vt(I);  % thermal voltage of the I^th Diode  element 
    Is = elementList.Diodes.Is(I);       % sturation current for the I^th Diode element
    
    % Votage at positive node of the diode :  
    if(n1~=0)
        Vp = X(n1);
    end
    % Votage at negative node of the diode :  
    if(n2~=0)
        Vn = X(n2);
    end
           
    if(n1~=0)&&(n2~=0)
        cur_dif = Is/Vt*(exp((Vp-Vn)/Vt));
        J(n1,n1) = J(n1,n1) + cur_dif;
        J(n1,n2) = J(n1,n2) - cur_dif;
        J(n2,n1) = J(n2,n1) - cur_dif;
        J(n2,n2) = J(n2,n2) + cur_dif;
    elseif(n1==0)
        cur_dif = Is/Vt*exp((-Vn)/Vt);
        J(n2,n2) = J(n2,n2) + cur_dif;
    elseif(n2==0)
        cur_dif = Is/Vt*exp(Vp/Vt);
        J(n1,n1) = J(n1,n1) + cur_dif;
    end  
end

% add BJts using this 
for I=1:elementList.BJTs.numElements
    
     
    %get Nodes Numbers
    cNode = elementList.BJTs.nodeNumbers(I,1);    % collector node.
    bNode = elementList.BJTs.nodeNumbers(I,2);    % base node.
    eNode = elementList.BJTs.nodeNumbers(I,3);    % emitter node
     
    % get other parameters of the BJTs 
    Vt = elementList.BJTs.Vt(I);  % thermal voltage of the I^th BJT
    Is = elementList.BJTs.Is(I);       % sturation current for the I^th BJT
    alphaR = elementList.BJTs.alphaR(I);
    alphaF =  elementList.BJTs.alphaF(I);

    Vc = X(cNode);
    Vb = X(bNode);
    Ve = X(eNode);
    
    J(cNode,cNode) = Is/Vt*exp((Vb-Vc)/Vt);
    J(cNode,bNode) = -Is/Vt*exp((Vb-Vc)/Vt)+alphaF*Is/Vt*exp((Vb-Ve)/Vt);
    J(cNode,eNode) = -alphaF*Is/Vt*exp((Vb-Ve)/Vt);

    J(bNode,cNode) = alphaR*Is/Vt*exp((Vb-Vc)/Vt)-Is/Vt*exp((Vb-Vc)/Vt);
    J(bNode,bNode) = Is/Vt*exp((Vb-Vc)/Vt)-alphaF*Is/Vt*exp((Vb-Ve)/Vt)+Is/Vt*exp((Vb-Ve)/Vt)-alphaR*Is/Vt*exp((Vb-Vc)/Vt);
    J(bNode,eNode) = alphaF*Is/Vt*exp((Vb-Ve)/Vt)-Is/Vt*exp((Vb-Ve)/Vt);

    J(eNode,cNode) = -alphaR*Is/Vt*exp((Vb-Vc)/Vt);
    J(eNode,bNode) = -Is/Vt*exp((Vb-Ve)/Vt)+alphaR*Is/Vt*exp((Vb-Vc)/Vt);
    J(eNode,eNode) = Is/Vt*exp((Vb-Ve)/Vt);

end 

