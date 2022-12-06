function diode(Name,Node1,Node2,varargin)
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
% set defaults

Is=1e-13;
Vt=0.025;


for I = 1:2:size(varargin,2)
    switch upper( varargin{I} )
        case {'IS'}
            Is = varargin{I+1};
        case {'VT'}
            Vt = varargin{I+1};  
        otherwise
            error([varargin{I} ' is an invalid optional argument  name in voltage source ', ...
                Name ]);
    end % switch
end



idx = elementList.Diodes.containerMap.Count + 1; % element Index
elementList.Diodes.numElements = idx;  %number of elements in Diodes

elementList.Diodes.Name{idx} = Name;  %name of the Diodes
elementList.Diodes.containerMap = addElement(elementList.Diodes.containerMap, Name); %store the Name in the container Map
elementList.Diodes.NodeNames{idx} = {Node1 , Node2}; %store the node names

[nodeMap, elementList.Diodes.nodeNumbers(idx,:)] = addNodes(elementList.Diodes.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers


elementList.Diodes.Is(idx) = Is;
elementList.Diodes.Vt(idx) = Vt;


% update the nodes count
elementList.numNodes = nodeMap.Count;

end
