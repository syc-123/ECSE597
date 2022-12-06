function gamma =  makeGamma(H)

Nh = 2*H+1; % number of fourier coefficients
gamma = zeros(Nh,Nh);

for I = 1:Nh
    gamma(I,1) = 1;
end

for I=1:H
    for J=1:Nh
        gamma(J,(2*I)) = cos(2*pi/Nh*I*(J-1));
        gamma(J,(2*I+1)) = sin(2*pi/Nh*I*(J-1));
    end

end