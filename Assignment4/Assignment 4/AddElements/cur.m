function cur(Name,Node1,Node2,varargin)
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

% defaults
VAL_DC = 0;
FREQUENCY = 0;
PHASE = 0;
AMPLITUDE = 0;
Nh = 0;
for I = 1:2:size(varargin,2)
    switch upper( varargin{I} )
        case {'TYPE'}
            TYPE = varargin{I+1};
        case {'VAL_DC'}
            VAL_DC = varargin{I+1}
        case {'AMPLITUDE'}
            AMPLITUDE = varargin{I+1};
        case {'FREQUENCY'}
            FREQUENCY = varargin{I+1};
        case {'NUMHARMONICS'}
            Nh = varargin{I+1};
        case {'PHASE'}
            PHASE = varargin{I+1} ;
        otherwise
            error([varargin{I} ' is an invalid optional argument  name in voltage source ', ...
                Name ]);
    end % switch
end


switch TYPE
    
    case 'DC'
        
        idx = elementList.DC_CurSources.containerMap.Count + 1; % element Index
        elementList.DC_CurSources.numElements = idx;  %number of elements in CurSources
        
        elementList.DC_CurSources.Name{idx} = Name;  %name of the inductor
        elementList.DC_CurSources.containerMap = addElement(elementList.DC_CurSources.containerMap, Name); %store the Name in the container Map
        elementList.DC_CurSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        [nodeMap, elementList.DC_CurSources.nodeNumbers(idx,:)] = addNodes(elementList.DC_CurSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.DC_CurSources.VAL_DC = VAL_DC;
        
        
    case  'AC'
        idx = elementList.AC_CurSources.containerMap.Count + 1; % element Index
        elementList.AC_CurSources.numElements = idx;  %number of elements in CurSources
        
        elementList.AC_CurSources.Name{idx} = Name;  %name of the inductor
        elementList.AC_CurSources.containerMap = addElement(elementList.AC_CurSources.containerMap, Name); %store the Name in the container Map
        elementList.AC_CurSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        
        [nodeMap, elementList.AC_CurSources.nodeNumbers(idx,:)] = addNodes(elementList.AC_CurSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.AC_CurSources.VAL_DC = VAL_DC;
        elementList.AC_CurSources.AMPLITUDE = AMPLITUDE;
        %          elementList.AC_CurSources.FREQUENCY = FREQUENCY;
        elementList.AC_CurSources.PHASE = PHASE;
        
        
    case 'TRANSIENT'
        
        idx = elementList.TR_CurSources.containerMap.Count + 1; % element Index
        elementList.TR_CurSources.numElements = idx;  %number of elements in CurSources
        
        elementList.TR_CurSources.Name{idx} = Name;  %name of the inductor
        elementList.TR_CurSources.containerMap = addElement(elementList.TR_CurSources.containerMap, Name); %store the Name in the container Map
        elementList.TR_CurSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        [nodeMap, elementList.TR_CurSources.nodeNumbers(idx,:)] = addNodes(elementList.TR_CurSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.TR_CurSources.VAL_DC = VAL_DC;
        elementList.TR_CurSources.AMPLITUDE = AMPLITUDE;
        elementList.TR_CurSources.FREQUENCY = FREQUENCY;
        elementList.TR_CurSources.PHASE = PHASE;
        
    case 'HB'
        
        idx = elementList.HB_CurSources.containerMap.Count + 1; % element Index
        elementList.HB_CurSources.numElements = idx;  %number of elements in CurSources
        
        elementList.HB_CurSources.Name{idx} = Name;  %name of the inductor
        elementList.HB_CurSources.containerMap = addElement(elementList.HB_CurSources.containerMap, Name); %store the Name in the container Map
        elementList.HB_CurSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        elementList.HB_CurSources.VAL_DC = VAL_DC;
        elementList.HB_CurSources.AMPLITUDE = AMPLITUDE;
        %             elementList.TR_CurSources.FREQUENCY = FREQUENCY;
        elementList.HB_CurSources.PHASE = PHASE;
        
        [nodeMap, elementList.HB_CurSources.nodeNumbers(idx,:)] = addNodes(elementList.HB_CurSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
    otherwise
        error([TYPE ' is an invalid "type" argument  name in voltage source ', ...
            Name ]);
end

% update the nodes count 
    elementList.numNodes = nodeMap.Count;

end
