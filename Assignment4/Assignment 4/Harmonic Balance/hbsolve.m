function [Xout] = hbsolve(Xguess,H)

%This function solves the Harmonic Balance Equations using Newton-Raphson
%method.

% Input: Xguess:  is the intial guess in the frequency domain (i.e.
% containing  the fourier coefficients)
%       H   :  number of Hamonics
global elementList

Nh = 2*H+1;
Gbar = makeHB_Gmat(H); % harmonic balance matrix Gbar
Cbar = makeHB_Cmat(H); % harmonic balance matrix Cbar
n = elementList.n; % size of regular MNA
nHB = Nh*n; % size of Harmonc Balance MNA


[~,~,Bbar]= makeBvector(H);    % ***---  Bbar is the SOURCE VECTOR FOR HARMONIC BALANCE---***

freq = elementList.HB_VolSources.FREQUENCY ;  % ***--- frequency of the Harmonic Balance source ---***