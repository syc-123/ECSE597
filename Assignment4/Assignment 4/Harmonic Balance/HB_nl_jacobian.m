function J = HB_nl_jacobian(Xb,H)

%  this function takes the time domain vector Xs as input and returns
%  Jacobian J in FREQUENCY domain.

% Inputs: 1.  Xs is the Newton-Raphson guess vector in time domain
%        2. H is the number of harmonics

% Output: J is the jacobain of the nonlinear vector in the frequency dommain.


global elementList

n = elementList.n;
Nh = 2*H+1; % number of fourier coefficients.

%% Fill in the J for Diodes
