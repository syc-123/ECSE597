function Cbar = makeHB_Cmat(H)
% Inputs: H is the number of harmonics
% Output: Cbar is the Harmonic balnce matrix.

global elementList
C= makeCmatrix;
Nh = 2*H+1; % number of fourier coefficients
Cbar = zeros(Nh*size(C,1));
for I = 1:size(C,1)
    for J = 1:size(C,1)
      
        cmul=1;
        cval = C(I,J);
        for k = 1:Nh
                        
            if rem(k,2) == 0
                Cbar(Nh*J-(Nh-k)+1,Nh*I-(Nh-k))= (Cbar(Nh*J-(Nh-k)+1,Nh*I-(Nh-k)) + (-1*cval*(cmul)));
                Cbar(Nh*J-(Nh-k),Nh*I-(Nh-k)+1)= (Cbar(Nh*J-(Nh-k),Nh*I-(Nh-k)+1) + (cval*(cmul)));
                cmul = cmul +1;
            end
        end
    end
end
end 