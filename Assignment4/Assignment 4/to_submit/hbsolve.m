function [Xout] = hbsolve(Xguess,H)

global elementList

maxerr = 1e-8;

Nh = 2*H+1;
Gbar = makeHB_Gmat(H); % harmonic balance matrix Gbar
Cbar = makeHB_Cmat(H); % harmonic balance matrix Cbar
n = elementList.n; % size of regular MNA
nHB = Nh*n; % size of Harmonc Balance MNA
gamma = makeGamma(H);

large_gamma = zeros(nHB,nHB);
for I=1:n
    large_gamma((I-1)*Nh+1:I*Nh,(I-1)*Nh+1:I*Nh) = gamma;
end

[~,~,Bbar]= makeBvector(H);   

while 1
    Fbar = HB_fvect((large_gamma*Xguess),H);
    J = HB_nl_jacobian((large_gamma*Xguess),H);
    phi_x = Gbar * Xguess + Cbar * Xguess + Fbar - Bbar;
    delta_x = -J\phi_x;
    Xguess = Xguess + delta_x;
    Xout = Xguess;
    if norm(delta_x) < maxerr
        break
    end

end


