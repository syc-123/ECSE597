function [D,S] = sens_differentiation_method(fpoints,eleNames,out)
% This function uses DIFFERENTIATION method to compute the sensitivity of the 
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

out_NodeNumber = getNodeNumber(outNode);

D = zeros(length(fpoints),length(eleNames));
S = zeros(length(fpoints),length(eleNames));

Gmat = makeGmatrix;
Cmat = makeCmatrix;
b = makeBvector;

%%
for I=1:length(eleNames)
    eleName = eleNames{I};
    A_d = zeros(size(Gmat,1),size(Gmat,2));

    switch upper(eleName(1)) 
        
        case 'R' 
            eleIdx = elementList.Resistors.containerMap(eleName);
            
            nodes = elementList.Resistors.nodeNumbers(eleIdx,:);

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

        x_d = A \ (-A_d * x);
    
        D(i,I) = x_d(out_NodeNumber,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S(i,I) = x_d(out_NodeNumber,1)/x(out_NodeNumber,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end

end

end