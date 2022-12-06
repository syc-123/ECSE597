% Testbench for Question 3 
 clear all


% load netlist 
Q3BEcircuit;
add_extraIndices;

t1 = 0; t2 = 15e-9; out = '3';
h1 = 1e-12; h2 = 30e-12; h3 = 60e-12;
% Find transient response using backward Euler using the three step sizes
[tpointsh1,tranrespBEh1] = BE_method(t2,h1,out);
[tpointsh2,tranrespBEh2] = BE_method(t2,h2,out);
[tpointsh3,tranrespBEh3] = BE_method(t2,h3,out);


% Find transient response using backward Euler using the three step sizes
[tpointsh1,tranrespFEh1] = FE_method(t2,h1,out);
[tpointsh2,tranrespFEh2] = FE_method(t2,h2,out);
[tpointsh3,tranrespFEh3] = FE_method(t2,h3,out);


close all
% Plot and compare results
figure(1)
hold off
clf
plot(1000*tpointsh1, 1000*tranrespBEh1,'b','LineWidth',2);
hold on
plot(1000*tpointsh2, 1000*tranrespBEh2,'r:','LineWidth',2);
plot(1000*tpointsh3, 1000*tranrespBEh3,'g:','LineWidth',2);
grid on
xlabel('Time (mSec)','FontSize',12)
ylabel('Output Voltage (mV)','FontSize',12)
legend({'Backward Euler h=1pSec','Backward Eurler h=30pSec', ...
    'Trapezoidal Rule h=60pSec'},'FontSize',12);



% Plot and compare results for Forward Euler
figure(2)
hold off
clf
plot(1000*tpointsh1, 1000*tranrespFEh1,'b','LineWidth',2);


plot(1000*tpointsh2, 1000*tranrespFEh2,'r:','LineWidth',2);
grid on
xlabel('Time (mSec)','FontSize',12)
ylabel('Output Voltage (mV)','FontSize',12)
legend({'Forward Euler h=1pSec','Forward Euler h=30pSec'},'FontSize',12);

% Plot forward Euler for h=60ps
figure(3)
hold off
clf
plot(1000*tpointsh3, 1000*tranrespFEh3,'g','LineWidth',2);
grid on
xlabel('Time (mSec)','FontSize',12)
ylabel('Output Voltage (mV)','FontSize',12)
legend({'Forward Euler h=60pSec'},'FontSize',12)
