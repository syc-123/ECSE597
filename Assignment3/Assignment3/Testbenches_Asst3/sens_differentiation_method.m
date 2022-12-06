function [D,S] = sens_differentiation_method(fpoints,eleNames,outNode)

global elementList

out_NodeNumber = getNodeNumber(outNode);

D = zeros(length(fpoints),length(eleNames));
S = zeros(length(fpoints),length(eleNames));

Gmat = makeGmatrix;
Cmat = makeCmatrix;
[bdc, b] = makeBvector;

%%
for I=1:length(eleNames)
    eleName = eleNames{I};
    G_d = zeros(size(Gmat));
    C_d = zeros(size(Gmat));

    switch upper(eleName(1)) 
        
        case 'R' 
            eleIdx = elementList.Resistors.containerMap(eleName);
            
            nodes = elementList.Resistors.nodeNumbers(eleIdx,:);

            val = 1/elementList.Resistors.value(eleIdx);

            if(nodes(1)~=0) && (nodes(2)~=0)
                G_d(nodes(1),nodes(1)) = 1;
                G_d(nodes(1),nodes(2)) = -1;
                G_d(nodes(2),nodes(1)) = -1;
                G_d(nodes(2),nodes(2)) = 1;
            elseif (nodes(1)==0) && (nodes(2)~=0)        
                G_d(nodes(2),nodes(2)) = 1;       
            elseif (nodes(1)~=0) && (nodes(2)==0)
                G_d(nodes(1),nodes(1)) = 1;
            end 

        case 'C' 

            eleIdx = elementList.Capacitors.containerMap(eleName);

            nodes = elementList.Capacitors.nodeNumbers(eleIdx,:);

            val = elementList.Capacitors.value(eleIdx);
            
            if(nodes(1)~=0) && (nodes(2)~=0)
                C_d(nodes(1),nodes(1)) = 1;
                C_d(nodes(1),nodes(2)) = -1;
                C_d(nodes(2),nodes(1)) = -1;
                C_d(nodes(2),nodes(2)) = 1;
            elseif (nodes(1)==0) && (nodes(2)~=0)
                C_d(nodes(2),nodes(2)) = 1;       
            elseif (nodes(1)~=0) && (nodes(2)==0)
                C_d(nodes(1),nodes(1)) = 1;
            end

    end


    for i = 1:length(fpoints)
        A = Gmat + 2 * pi * fpoints(i) * 1i * Cmat;
        A_d = G_d + 2 * pi * fpoints(i) * 1i * C_d;
        x = A \ b;

        x_d = A \ (-A_d * x);
    
        D(i,I) = x_d(out_NodeNumber,1);
        S(i,I) = val/x(out_NodeNumber,1)*D(i,I);
    
    end

end

end