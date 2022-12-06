% Testbench for Q2

% load netlist for filter with time domain source
Circuit_chebychev_filter_TD;
add_extraIndices;

% Find transient response using backward Euler step size of 1 microSec
[tpointsh1,tranrespBEh1] = BE_method(20e-3,1e-6,out);

% Find transient response using BE with step size of 10 microSec
[tpointsh2,tranrespBEh2] = BE_method(20e-3,1e-5,out);

% Find transient response using Trapezoidal rule with step size of 10 microsec
[tpointsh2,tranrespTRh2] = Trapezoidal_method(20e-3,1e-5,out);



% Plot and compare results
figure(1)
hold off
clf
plot(1000*tpointsh1, 1000*tranrespBEh1,'b','LineWidth',2);
hold on
plot(1000*tpointsh2, 1000*tranrespBEh2,'r:','LineWidth',2);
plot(1000*tpointsh2, 1000*tranrespTRh2,'g:','LineWidth',2);
grid on
xlabel('Time (mSec)')
ylabel('Output Voltage (mV)')
legend({'Backward Euler h=1\muSec','Backward Eurler h=10\muSec', ...
    'Trapezoidal Rule h=10\muSec'});
