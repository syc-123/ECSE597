function [D,S] = sens_perturbation_method(fpoints,eleNames,outNode)
% This function uses DIFFERENCE method to compute the sensitivity of the 
% output node with respect to all the parameters.
%Inputs: 1. fpoints: contains the frequency points at which the sensitivity is
%                 is required.
%        2. eleNames: is a cell array contains the names of the elements.
%        3. outNode: is the node for which the sensitivity is required.

%Output:  1.D:  is a matrix.  It should contain the ABSOLUTE sensitivity of the
%              outNode at all fpoints for all elements.
%              One way to fill store senstivity in  D is to  add sensitivity of 
%              a given element in one column of D for all fpoints.
%              In this case if there are F number of frequency points and 
%              P number of elements in eleNames, then the size of matrix D
%              will be FxP.
%
%         2. S: is  matrix. It should contain the RELATIVE sensitivity of 
%               outNode for all the elements in eleNames. It can be filled
%               similar to matrix D.
%
%
%Note: There are multiple ways to add the code for computing  sensitivity..
%      The function stub provided is just an example. You can modify the 
%      in function in any fashion. Even you can split this function in to
%      multiple functions.
%--------------------------------------------------------------------------
global elementList

delta = 1e-7;

out_NodeNumber = getNodeNumber(outNode);

D = zeros(length(fpoints),length(eleNames));
S = zeros(length(fpoints),length(eleNames));

Gmat = makeGmatrix;
Cmat = makeCmatrix;
b = makeBvector;
f_response = zeros(1,length(fpoints));

for i = 1:length(fpoints)
    
    A = Gmat + 2 * pi * fpoints(i) * 1i * Cmat;
    
    x = A \ b;
    
    f_response(1,i) = x(out_NodeNumber,1);
    
end



%%
for I=1:length(eleNames)
    eleName = eleNames{I};
    G_d = zeros(size(Gmat,1),size(Gmat,2));
    C_d = zeros(size(Cmat,1),size(Cmat,2));
    f_response_new = zeros(1,length(fpoints));

    switch upper(eleName(1))
        
        case 'R'

            eleIdx = elementList.Resistors.containerMap(eleName);

            nodes = elementList.Resistors.nodeNumbers(eleIdx,:);

            val = elementList.Resistors.value(eleIdx);

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
        
        A_new = (Gmat + G_d) + 2 * pi * fpoints(i) * 1i * (Cmat + C_d);

        x_new = A_new \ b;
    
        f_response_new(1,i) = x_new(out_NodeNumber,1);

        D(i,I) = f_response(1,i) - f_response_new(1,i);
        S(i,I) = D(i,I)/f_response(1,i);
    
    end

end

end