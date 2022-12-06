function Btr = makeBt(time)
% This function adds the sinusoidal source for Bvector for transient
% analysis at a given value of time.
% Output : B is the ouptut time vector.
% Input: time is the time in seconds at which the sinusoidal source is to
% be evaluated.
%
% Note: This function only adds the  Bvector stamp for Transient sources.
%       The contribution of the Transient voltage sourcees needs to be
%       added in the G matrix. Also, if the DC sources are present in the
%       circuit netlist, then their contribution in Bvector need to be added
%       using the function named makeBvector.m. The function makeBvector.m
%       is provided along with this Assignment.
%
%--------------------------------------------------------------------------

%%
global elementList;
Btr = zeros(elementList.n,1); % declare the transient B vector.


for I = 1:elementList.TR_VolSources.numElements
    % get the curIdx. The name for current index may be different for you,
    %  plese check your code.
    
    curIdx = elementList.TR_VolSources.curIdx(I);
    
    Amp = elementList.TR_VolSources.AMPLITUDE(I);
    Phase = elementList.TR_VolSources.PHASE(I);
    Val_DC = elementList.TR_VolSources.VAL_DC(I);
    Freq = elementList.TR_VolSources.FREQUENCY(I);
    
    %     add it for type sinusoidal sources
    Btr(curIdx,1) = Val_DC + Amp*sin(2*pi*Freq*time);
    
end


for I = 1:elementList.TR_CurSources.numElements
    % get the curIdx. The name for current index may be different for you,
    %  plese check your code.
    
    %     curIdx = elementList.TR_VolSources.curIdx(I);
    %
    Amp = elementList.TR_CurSources.AMPLITUDE(I);
    Val_DC = elementList.TR_CurSources.VAL_DC(I);
    Freq = elementList.TR_CurSources.FREQUENCY(I);
    %
    %     %     add it for type sinusoidal sources
    %     Btr(curIdx,1) = Val_DC + Amp*sin(2*pi*Freq*time);
    
    nodes = elementList.TR_CurSources.nodeNumbers(I,:);
    val = Val_DC + Amp*sin(2*pi*Freq*time);
    
    if(nodes(1)~=0)
        Btr(nodes(1),1) = val;
    end
    
    if(nodes(2)~=0)
        Btr(nodes(2),1) = -val;
    end
    
end