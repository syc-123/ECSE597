clear all


global elementList

Half_wave_rectifier_TR; %circuit netlist
%% Harmonic Balance Simulation
% CHOOSE THE SUITABLE number of harmonics
          
% H =12                                         % H is the number of harmonics
% Nh = 2*H+1;                                        % Nh is the number of fourier coefficients.
% Xguess = zeros(elementList.n*Nh ,1);                     % Make Xguess.

% [Xout] = hbsolve(Xguess,H)                         % Solve the Harmonic Balance system.

[tpoints, Xout]= nl_BE_method(0.05,1e-4);
%% plot here
close all
hold on 
plot(tpoints,Xout(1,:),'b-','LineWidth',2,'displayname','Voltage  at Node 1')
plot(tpoints,Xout(2,:),'k-','LineWidth',2,'displayname','Voltage  at Node 2')
plot(tpoints,Xout(3,:),'r-','LineWidth',2,'displayname','Voltage  at Node 3')
grid on 
ylabel('Voltage [V]')
xlabel('Time [s]')
title( 'Node voltages of Half-wave rectifier using Transient Analysis')
% legend('Voltage  at Node 1','Voltage  at Node 2' , 'Voltage  at Node 3' )
legend
% yticks([-1.5:0.2:1.5])