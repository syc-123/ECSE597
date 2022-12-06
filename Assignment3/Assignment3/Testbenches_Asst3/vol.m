function vol(Name,Node1,Node2,varargin)
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
            VAL_DC = varargin{I+1};
        case {'AMPLITUDE'}
            AMPLITUDE = varargin{I+1};
        case {'FREQUENCY'}
            FREQUENCY = varargin{I+1};
        case {'NUMHARMONICS'}
            Nh = varargin{I+1};
        case {'PHASE'}
            PHASE = varargin{I+1};
        otherwise
            error([varargin{I} ' is an invalid optional argument  name in voltage source ', ...
                Name ]);
    end % switch
end


switch upper(TYPE)
    
    case 'DC'
        
        idx = elementList.DC_VolSources.containerMap.Count + 1; % element Index
        elementList.DC_VolSources.numElements = idx;  %number of elements in VolSources
        
        elementList.DC_VolSources.Name{idx} = Name;  %name of the inductor
        elementList.DC_VolSources.containerMap = addElement(elementList.DC_VolSources.containerMap, Name); %store the Name in the container Map
        elementList.DC_VolSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        [nodeMap, elementList.DC_VolSources.nodeNumbers(idx,:)] = addNodes(elementList.DC_VolSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.DC_VolSources.VAL_DC(idx) = VAL_DC;
    case  'AC'
        idx = elementList.AC_VolSources.containerMap.Count + 1; % element Index
        elementList.AC_VolSources.numElements = idx;  %number of elements in VolSources
        
        elementList.AC_VolSources.Name{idx} = Name;  %name of the inductor
        elementList.AC_VolSources.containerMap = addElement(elementList.AC_VolSources.containerMap, Name); %store the Name in the container Map
        elementList.AC_VolSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        [nodeMap, elementList.AC_VolSources.nodeNumbers(idx,:)] = addNodes(elementList.AC_VolSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.AC_VolSources.VAL_DC(idx) = VAL_DC;
        elementList.AC_VolSources.AMPLITUDE(idx)= AMPLITUDE;
        %          elementList.AC_VolSources.FREQUENCY = FREQUENCY;
        elementList.AC_VolSources.PHASE = PHASE;
        
    case 'TRANSIENT'
        
        idx = elementList.TR_VolSources.containerMap.Count + 1; % element Index
        elementList.TR_VolSources.numElements = idx;  %number of elements in VolSources
        
        elementList.TR_VolSources.Name{idx} = Name;  %name of the inductor
        elementList.TR_VolSources.containerMap = addElement(elementList.TR_VolSources.containerMap, Name); %store the Name in the container Map
        elementList.TR_VolSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        
        [nodeMap, elementList.TR_VolSources.nodeNumbers(idx,:)] = addNodes(elementList.TR_VolSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.TR_VolSources.VAL_DC = VAL_DC;
        elementList.TR_VolSources.AMPLITUDE = AMPLITUDE;
        elementList.TR_VolSources.FREQUENCY = FREQUENCY;
        elementList.TR_VolSources.PHASE = PHASE;
        
    case 'HB'
        
        idx = elementList.HB_VolSources.containerMap.Count + 1; % element Index
        elementList.HB_VolSources.numElements = idx;  %number of elements in VolSources
        
        elementList.HB_VolSources.Name{idx} = Name;  %name of the inductor
        elementList.HB_VolSources.containerMap = addElement(elementList.HB_VolSources.containerMap, Name); %store the Name in the container Map
        elementList.HB_VolSources.NodeNames{idx} = {Node1 , Node2}; %store the node names
        [nodeMap, elementList.HB_VolSources.nodeNumbers(idx,:)] = addNodes(elementList.HB_VolSources.NodeNames{idx},nodeMap); % add nodes to Map and get node numbers
        
        elementList.HB_VolSources.VAL_DC = VAL_DC;
        elementList.HB_VolSources.AMPLITUDE = AMPLITUDE;
        %             elementList.HB_VolSources.FREQUENCY = FREQUENCY;
        elementList.HB_VolSources.PHASE = PHASE;
        
end

% update the nodes count 
    elementList.numNodes = nodeMap.Count;

end
