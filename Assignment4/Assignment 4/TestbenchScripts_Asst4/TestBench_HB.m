clear all


global elementList

Half_wave_rectifier_HB; %circuit netlist
%% Harmonic Balance Simulation


% CHOOSE THE SUITABLE number of harmonics!!!

H =1                                         % H is the number of harmonics
Nh = 2*H+1;                                  % Nh is the number of fourier coefficients.
Xguess = zeros(elementList.n*Nh ,1);         % Make Xguess.

[Xout] = hbsolve(Xguess,H);                   % Solve the Harmonic Balance system.


%% plot the input and output node results(already Done)
Gamma = makeGamma(H);
outNode =3;                                                  % output node
XoutNode = Xout( (outNode-1)*Nh +1: outNode*Nh);             % get the fourier coefficients at output node
XoutNode = Gamma*XoutNode;                                   % Convert the fourier coefficeints to time domain .

InputNode = 1;                                              % input node
XinNode = Xout( (InputNode-1)*Nh +1: InputNode*Nh);         % get the fourier coefficients at input node
XinNode = Gamma*XinNode;                                    % Convert the fourier coefficeints to time domain .

XinNode = [XinNode; XinNode];
XoutNode =[XoutNode;XoutNode];

node2 = 2;
Xnode2 = Xout( (node2-1)*Nh +1: node2*Nh);             % get the fourier coefficients at output node

Xnode2= Gamma*Xnode2;                                   % Convert the fourier coefficeints to time domain .
Xnode2 = [Xnode2;Xnode2];

time = 0:1/(60*Nh):1/30;
time = time(1:end-1);
% plot here
close all
hold on
plot(time,XinNode,'b-','LineWidth',2)
plot(time,Xnode2,'k-','LineWidth',2)
plot(time,XoutNode,'r-','LineWidth',2)
grid on
ylabel('Voltage [V]')
xlabel('Time [s]')
title( 'Voltage at input and output of Half-wave rectifier')
legend('Voltage  at Node 1','Voltage  at Node 2' , 'Voltage  at Node 3' )

yticks([-1.5:0.2:1.5])