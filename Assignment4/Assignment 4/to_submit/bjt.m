function bjt(Name,Node1,Node2,Node3,varargin)
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
TYPE = 'NPN';
Is=1e-13;
Vt=0.025;
alphaR = 0.1;
alphaF = 0.99;

for I = 1:2:size(varargin,2)
    switch upper( varargin{I} )
        case {'TYPE'}
            TYPE = varargin{I+1};
        case {'IS'}
            Is = varargin{I+1};
        case {'VT'}
            Vt = varargin{I+1};
        case {'ALPHAF'}
            alphaF = varargin{I+1};
        case {'ALPHAR'}
            alphaR = varargin{I+1};
        
        otherwise
            error([varargin{I} ' is an invalid optional argument  name in voltage source ', ...
                Name ]);
    end % switch
end



idx = elementList.BJTs.containerMap.Count + 1; % element Index
elementList.BJTs.numElements = idx;  %number of elements in BJTs

elementList.BJTs.Name{idx} = Name;  %name of the BJT
elementList.BJTs.containerMap = addElement(elementList.BJTs.containerMap, Name); %store the Name in the container Map
elementList.BJTs.NodeNames{idx} = {Node1 , Node2,Node3}; %store the node names

[nodeMap, elementList.BJTs.nodeNumbers(idx,:)] = addNodes(elementList.BJTs.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers


elementList.BJTs.TYPE{idx} =  TYPE;
elementList.BJTs.Is(idx) = Is;
elementList.BJTs.Vt(idx) = Vt;
elementList.BJTs.alphaR(idx) = alphaR;
elementList.BJTs.alphaF(idx) = alphaF;

% update the nodes count
elementList.numNodes = nodeMap.Count;

end
