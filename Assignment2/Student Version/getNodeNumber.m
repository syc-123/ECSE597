function nodeNumbers = getNodeNumber(nodeNames)
% Author: Karan 
% Date: 13/09/2022

% Input: nodeNames is a cell array of nodeNames 
% Ouput: nodeNumber is a vector of associates node numbers 

% Example 1: for a single  node named 'n1'  call this function as 
%             nodeNumbers = getNodeNumber('n1')




global nodeMap

% nodeNumbers = zeros(1,size(nodeNames,2));
% for I=1:size(nodeNames,2)
    
    if(strcmp(nodeNames,'0')  || strcmp(nodeNames,'gnd') )
        nodeNumbers = 0;
    elseif( nodeMap.isKey(nodeNames))
        nodeNumbers = nodeMap(nodeNames);
    else
         error(['The node name ' nodeNames  ' is not assigned in the netlist.'...
             'Check your input.']);

    end
    
% end