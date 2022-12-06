function res(Name,Node1,Node2,value)
% Author: Karan 
% Date: 13/09/2022
    % Name:  type: Char
    %        Name of the element, this must be unique for each element.
    
    % Node1, Node2:  type: Char
    %        Name of the node1 and  node2.
    
    % Value: element value ( inductance, capacitance, resistance etc.)  
                % for  inductors unit of value is Henr
                % for  resistors unit of value is Ohm
                % for  capacitors unit of value is Farad
                % for  voltage source unit of value is Volts 
                % for current sources unit of value is Ampere, etc.
                
    global elementList  nodeMap

    idx = elementList.Resistors.containerMap.Count + 1; % element Index
    elementList.Resistors.numElements = idx;  %number of elements in Resistors
    
    elementList.Resistors.Name{idx} = Name;  %name of the inductor
    elementList.Resistors.containerMap = addElement(elementList.Resistors.containerMap, Name); %store the Name in the container Map 
    elementList.Resistors.NodeNames{idx} = {Node1 , Node2}; %store the node names  

    
    [nodeMap, elementList.Resistors.nodeNumbers(idx,:)] = addNodes(elementList.Resistors.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers 
    
    % store the value 
    elementList.Resistors.value(idx) =  value;  %store the value 
    
    % update the nodes count 
    elementList.numNodes = nodeMap.Count;
end 

