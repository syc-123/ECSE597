function [D,S] = adjoint_sensitivity(fpoints,eleNames,outNode)

global elementList

%%
out_NodeNumber = getNodeNumber(outNode);

D = zeros(length(fpoints),length(eleNames));
S = zeros(length(fpoints),length(eleNames));

Gmat = makeGmatrix;
Cmat = makeCmatrix;
[bdc, b] = makeBvector;

d = zeros(size(Gmat,2),1);
d(out_NodeNumber,1) = 1;

%%
for I=1:length(eleNames)
    eleName = eleNames{I};
    A_d = zeros(size(Gmat));

    switch upper(eleName(1)) 
        
        case 'R' 
            eleIdx = elementList.Resistors.containerMap(eleName);
            
            nodes = elementList.Resistors.nodeNumbers(eleIdx,:);

            val = elementList.Resistors.value(eleIdx);

            if(nodes(1)~=0) && (nodes(2)~=0)
                A_d(nodes(1),nodes(1)) = 1;
                A_d(nodes(1),nodes(2)) = -1;
                A_d(nodes(2),nodes(1)) = -1;
                A_d(nodes(2),nodes(2)) = 1;
            elseif (nodes(1)==0) && (nodes(2)~=0)        
                A_d(nodes(2),nodes(2)) = 1;       
            elseif (nodes(1)~=0) && (nodes(2)==0)
                A_d(nodes(1),nodes(1)) = 1;
            end 

        case 'C' 

            eleIdx = elementList.Capacitors.containerMap(eleName);

            nodes = elementList.Capacitors.nodeNumbers(eleIdx,:);

            val = elementList.Capacitors.value(eleIdx);
            
            if(nodes(1)~=0) && (nodes(2)~=0)
                A_d(nodes(1),nodes(1)) = 1;
                A_d(nodes(1),nodes(2)) = -1;
                A_d(nodes(2),nodes(1)) = -1;
                A_d(nodes(2),nodes(2)) = 1;
            elseif (nodes(1)==0) && (nodes(2)~=0)
                A_d(nodes(2),nodes(2)) = 1;       
            elseif (nodes(1)~=0) && (nodes(2)==0)
                A_d(nodes(1),nodes(1)) = 1;
            end

    end


    for i = 1:length(fpoints)
        A = Gmat + 2 * pi * fpoints(i) * 1i * Cmat;
        
        x = A \ b;
        x_a = transpose(A) \ (-d);
        

    
        D(i,I) = transpose(x_a) * A_d * x;
        S(i,I) = val/ x(out_NodeNumber,1)*D(i,I);
        
    
    end

end

end