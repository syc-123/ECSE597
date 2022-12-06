function Fbar = HB_fvect(Xs,H)
%Xs is in time-domain
global elementList

n = elementList.n;
Nh = 2*H+1;
nHB = Nh*n;

Is=1e-13;
Vt=0.025;

gamma = makeGamma(H);
F = zeros(nHB,1);

for K=1:elementList.Diodes.numElements
    nodes = elementList.Diodes.nodeNumbers(K,:);
   
    n1 = nodes(K,1);
    n2 = nodes(K,2);
     
    if(nodes(1)~=0 && nodes(2)~=0)
        for I=1:Nh
            Vd = Xs((n1-1)*Nh+I,1)-Xs((n2-1)*Nh+I,1);
            F((n1-1)*Nh+I,1) = F((n1-1)*Nh+I,1) + Is*(exp(Vd/Vt)-1);
            F((n2-1)*Nh+I,1) = F((n2-1)*Nh+I,1) - Is*(exp(Vd/Vt)-1);
        end 
    elseif (nodes(1)==0 && nodes(2)~=0)
        for I=1:Nh
            Vd = -Xs((n2-1)*Nh+I,1);
            F((n2-1)*Nh+I,1) = F((n2-1)*Nh+I,1) - Is*(exp(Vd/Vt)-1);
        end    
    elseif(nodes(1)~=0 && nodes(2)==0)
        for I=1:Nh
            Vd = Xs((n1-1)*Nh+I,1);
            F((n1-1)*Nh+I,1) = F((n1-1)*Nh+I,1) + Is*(exp(Vd/Vt)-1);
        end 
    end
end

large_gamma = zeros(nHB,nHB);
for I=1:n
    large_gamma((I-1)*Nh+1:I*Nh,(I-1)*Nh+1:I*Nh) = gamma;
end

Fbar = large_gamma\F;