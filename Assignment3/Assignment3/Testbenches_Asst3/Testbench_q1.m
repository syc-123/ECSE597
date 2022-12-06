% Testbench Q1
clear all
Circuit_chebychev_filter_TD;   % load netlist for filter with time domain source
add_extraIndices

% Find transient response using backward Euler
tEnd =20e-3;
h = 1e-6;
out = '10';

% compute the transient reposnse using BE
[tpoints, tranresp]= BE_method(tEnd,h, out);  % feel free to modify the function



Circuit_chebychev_filter_freq % Load chebychev filter with AC input
add_extraIndices

% Set frequency point at 2000Hz
fpoint = 2000;

% call the fsovle function to compute the frequency response at fpoint
Fresp = fsolve(fpoint , out);   % modify it according to your fsolve. 


%% Sample the steady state output in the time domain
steadystate = abs(Fresp)*sin(2*pi*fpoint*tpoints + angle(Fresp));

%plots
figure(1)
hold off
clf
plot(1000*tpoints, 1000*tranresp,'b','LineWidth',2);
hold on
plot(1000*tpoints, 1000*steadystate,'r:','LineWidth',2);
grid on
xlabel('Time (mSec)')
ylabel('Output Voltage (mV)')
legend('Transient', 'Steady State');




