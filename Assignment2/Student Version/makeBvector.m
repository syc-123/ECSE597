function [Bdc, Bac] = makeBvector()
% This function adds the contribution of input sources in the Bvector

% Bdc is the MNA RHS vector containing the contribution of  DC sources 
% Bac is the MNA RHS vector containing the contribution of  AC sources 
global elementList 

 
Bdc = sparse(elementList.n,1);
Bac = sparse(elementList.n,1);

% add stamps here...
for I=1:elementList.DC_VolSources.numElements
    % acess the current  index 
    nX = elementList.DC_VolSources.currIndex(I);
    Bdc(nX,1) =  elementList.DC_VolSources.VAL_DC(I);
end 

% add stamps here...
for I=1:elementList.AC_VolSources.numElements
    % acess the current  index 
    nX = elementList.AC_VolSources.curIdx(I);
    Bac(nX,1) =  elementList.AC_VolSources.AMPLITUDE(I);
end 

