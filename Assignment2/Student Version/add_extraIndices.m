function add_extraIndices()
%Author:Karan
%Date: 22/09/2022
%in this function we will add extra idices which are required by certain
% circuit elements 

global elementList
elementList.n = elementList.numNodes; % MNA size

% add extra indices for DC voltages
for I=1:elementList.DC_VolSources.numElements
    % asign extra index
    % for this we add a new field named curIdx in DC_VolSources
    
    elementList.DC_VolSources.curIdx(I) = elementList.n+1;
    
    % update  MNA size
    elementList.n =elementList.n +1;
end 

% add extra indices for AC voltages
% add extra indices for DC voltages
for I=1:elementList.AC_VolSources.numElements
    % asign extra index
    % for this we add a new field named curIdx in DC_VolSources
    
    elementList.AC_VolSources.curIdx(I) = elementList.n+1;
    
    % update  MNA size
    elementList.n =elementList.n +1;
end 
% add extra indices for Transient voltages

% add extra indices for Inductors 
