clear all
clc 


sallenkey_lp; % netlist.

add_extraIndices;

Gmat = makeGmatrix;
Cmat = makeCmatrix;

out = 'n4';
fpoints = linspace(0,1000,100);
r = fsolve(fpoints, out);

close all
figure(1)
semilogx(fpoints,20*log10(abs(r)),'linewidth',2,'displayname', ['\mid V_{' out '} \mid'])
ylabel(['\mid V_{' out '} \mid (dB)'])
xlabel('Frequency (Hz)')
legend
grid on

%%
elementNames = {'R3','R4','R1','C1','R2','C2'};
[Ddelta,Sdelta] = sens_perturbation_method(fpoints, elementNames,out);
[Ddiff,Sdiff] = sens_differentiation_method(fpoints,elementNames,out);
[Dadj,Sadj] = sens_differentiation_method(fpoints,elementNames,out);

for I= 1:length(elementNames)
figure(I+1)
semilogx(fpoints,abs(Sdelta(:,I)),'b-','linewidth',2,'displayname', [ 'Relative Sens. pertubation  for '  elementNames{I}])
hold on
semilogx(fpoints,abs(Sdiff(:,I)),'r--','linewidth',2,'displayname', ['Relative Sens. differentiation  for '  elementNames{I}])
semilogx(fpoints,abs(Sadj(:,I)),'g:','linewidth',1.5,'displayname', ['Relative Sens. adjoint  for '  elementNames{I}])
xlabel('Frequency (Hz)')
legend
grid on
end 

elementNames1 = {'R1','C1'};

for I= 1:length(elementNames1)
figure(I+7)
semilogx(fpoints,abs(Ddelta(:,I)),'b-','linewidth',2,'displayname', [ 'Abs Sens. pertubation  for '  elementNames1{I}])
hold on
semilogx(fpoints,abs(Ddiff(:,I)),'r--','linewidth',2,'displayname', ['Abs Sens. differentiation  for '  elementNames1{I}])
semilogx(fpoints,abs(Dadj(:,I)),'g:','linewidth',1.5,'displayname', ['Abs Sens. adjoint  for '  elementNames1{I}])
xlabel('Frequency (Hz)')
legend
grid on
end 