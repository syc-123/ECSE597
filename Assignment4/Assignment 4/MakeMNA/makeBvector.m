function [Bdc, Bac, Bhb] = makeBvector(H)
% This function adds the contribution of input sources in the Bvector
% Input:
 % H  is the number of Harmonics. H is only needed for the case of HB
 % sources.

% Outputs:
% Bdc is the MNA RHS vector containing the contribution of  DC sources 
% Bac is the MNA RHS vector containing the contribution of  AC sources 
% Bhb is the MNA RHS vector containing the contribution of  HB sources 

global elementList 

 
Bdc = sparse(elementList.n,1);
Bac = sparse(elementList.n,1);


% add stamps for DC sources here...
for I=1:elementList.DC_VolSources.numElements
    % acess the current  index 
    nX = elementList.DC_VolSources.curIdx(I);
    Bdc(nX,1) =  elementList.DC_VolSources.VAL_DC(I);
end 

% add stamps for AC here...
for I=1:elementList.AC_VolSources.numElements
    % acess the current  index 
    nX = elementList.AC_VolSources.curIdx(I);
    Bac(nX,1) =  elementList.AC_VolSources.AMPLITUDE(I);
end 

% add stamps for HB here...
if(nargin==1)
    % if there is a argument. then it is Harmonic Balance source is present
    Nh = 2*H+1;
    Bhb = zeros(elementList.n*Nh,1);
    
    for I=1:elementList.HB_VolSources.numElements
        % acess the current  index
        nX = elementList.HB_VolSources.curIdx(I);
        
        loc_sin =2 ;%location of sine
        %  this is added for the sinusoidal source.
        Bhb((nX-1)*Nh+1 + loc_sin,1) =  elementList.HB_VolSources.AMPLITUDE(I);
    end
end


