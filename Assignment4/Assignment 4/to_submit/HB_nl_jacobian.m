function J = HB_nl_jacobian(Xs,H)

%Xs is in time-domain

global elementList

n = elementList.n;
Nh = 2*H+1; % number of fourier coefficients.
nHB = Nh*n;

Is=1e-13;
Vt=0.025;

Gbar = makeHB_Gmat(H);
Cbar = makeHB_Cmat(H);
gamma = makeGamma(H);

large_gamma = zeros(nHB,nHB);
for I=1:n
    large_gamma((I-1)*Nh+1:I*Nh,(I-1)*Nh+1:I*Nh) = gamma;
end

%% Fill in the J for Diodes
J_t = zeros(nHB,nHB);

for I=1:elementList.Diodes.numElements
    nodes = elementList.Diodes.nodeNumbers(I,:);
   
    n1 = nodes(I,1);
    n2 = nodes(I,2);
     
    if(nodes(1)~=0 && nodes(2)~=0)
        for I=1:Nh
            Vd = Xs((n1-1)*Nh+I,1)-Xs((n2-1)*Nh+I,1);
            J_t((n1-1)*Nh+I,(n1-1)*Nh+I) = J_t((n1-1)*Nh+I,(n1-1)*Nh+I) + Is/Vt*exp(Vd/Vt);
            J_t((n2-1)*Nh+I,(n2-1)*Nh+I) = J_t((n2-1)*Nh+I,(n2-1)*Nh+I) + Is/Vt*exp(Vd/Vt);
        end 
    elseif (nodes(1)==0 && nodes(2)~=0)
        for I=1:Nh
            Vd = -Xs((n2-1)*Nh+I,1);
            J_t((n2-1)*Nh+I,(n2-1)*Nh+I) = J_t((n2-1)*Nh+I,(n2-1)*Nh+I) + Is/Vt*exp(Vd/Vt);
        end    
    elseif(nodes(1)~=0 && nodes(2)==0)
        for I=1:Nh
            Vd = Xs((n1-1)*Nh+I,1);
            J_t((n1-1)*Nh+I,(n1-1)*Nh+I) = J_t((n1-1)*Nh+I,(n1-1)*Nh+I) + Is/Vt*exp(Vd/Vt);
        end 
    end
end

J = Gbar + Cbar + large_gamma\J_t*large_gamma;