function Gbar = makeHB_Gmat(H)
% this function uses the MNA matrix, G, to compute the HB matrix Gbar 

% input: H is the number of harmonics.
% output: Gbar is the  Harnmonic balance matrix G,

global elementList 
G= makeGmatrix;
Nh = 2*H+1; % number of fourier coefficients 
Gbar = zeros(elementList.n*Nh);
for I=1:size(G,1)
     for J =1:size(G,2)

          if(G(I,J)~=0)
            Gbar(Nh*(I-1) + 1 : I*Nh,  Nh*(J-1) + 1 : J*Nh) = G(I,J)*eye(Nh);
          end 
     end 
end 