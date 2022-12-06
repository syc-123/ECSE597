function [D,S] = sens_perturbation_method(fpoints,eleNames,outNode)

global elementList

delta = 1e-7;

out_NodeNumber = getNodeNumber(outNode);

D = zeros(length(fpoints),length(eleNames));
S = zeros(length(fpoints),length(eleNames));

Gmat = makeGmatrix;
Cmat = makeCmatrix;
[bdc, b] = makeBvector;

%%
for I=1:length(eleNames)
    eleName = eleNames{I};
    G_d = sparse(size(Gmat,1),size(Gmat,2));
    C_d = sparse(size(Cmat,1),size(Cmat,2));
    
    switch upper(eleName(1))
        
        case 'R'

            eleIdx = elementList.Resistors.containerMap(eleName);
            
            nodes = elementList.Resistors.nodeNumbers(eleIdx,:);

            val = 1/elementList.Resistors.value(eleIdx);

            if(nodes(1)~=0) && (nodes(2)~=0)
                G_d(nodes(1),nodes(1)) = G_d(nodes(1),nodes(1)) + delta*val;
                G_d(nodes(1),nodes(2)) = G_d(nodes(1),nodes(2)) - delta*val;
                G_d(nodes(2),nodes(1)) = G_d(nodes(2),nodes(1)) - delta*val;
                G_d(nodes(2),nodes(2)) = G_d(nodes(2),nodes(2)) + delta*val;
            elseif (nodes(1)==0) && (nodes(2)~=0)        
                G_d(nodes(2),nodes(2)) = G_d(nodes(2),nodes(2)) + delta*val;       
            elseif (nodes(1)~=0) && (nodes(2)==0)
                G_d(nodes(1),nodes(1)) = G_d(nodes(1),nodes(1)) + delta*val;
            end 

        case 'C'

            eleIdx = elementList.Capacitors.containerMap(eleName);

            nodes = elementList.Capacitors.nodeNumbers(eleIdx,:);

            val = elementList.Capacitors.value(eleIdx);
 
            if(nodes(1)~=0) && (nodes(2)~=0)
                C_d(nodes(1),nodes(1)) = C_d(nodes(1),nodes(1)) + delta*val;
                C_d(nodes(1),nodes(2)) = C_d(nodes(1),nodes(2)) - delta*val;
                C_d(nodes(2),nodes(1)) = C_d(nodes(2),nodes(1)) - delta*val;
                C_d(nodes(2),nodes(2)) = C_d(nodes(2),nodes(2)) + delta*val;
            elseif (nodes(1)==0) && (nodes(2)~=0)
                C_d(nodes(2),nodes(2)) = C_d(nodes(2),nodes(2)) + delta*val;       
            elseif (nodes(1)~=0) && (nodes(2)==0)
                C_d(nodes(1),nodes(1)) = C_d(nodes(1),nodes(1)) + delta*val;
            end

    end

    
    for i = 1:length(fpoints)
        A = Gmat + 2 * pi * fpoints(i) * 1i * Cmat;
        A_new = (Gmat + G_d) + 2 * pi * fpoints(i) * 1i * (Cmat + C_d);
        x = A \ b;
        x_new = A_new \ b;
    
        D(i,I) = (x_new(out_NodeNumber,1)-x(out_NodeNumber,1))/(delta*val);
        S(i,I) = val/x(out_NodeNumber,1)*D(i,I);
    
    end

end

end