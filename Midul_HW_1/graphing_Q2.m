function [t,Vm] = graphing_Q2(Er,I0,T,R,Cm,tau)
% This function computes Tm in accordance to the solution to Q2a
% Define needed parameters and t vectors:
A = 4*pi*(R^2); w = 2/T; c=(I0/(A*Cm))*((pi*w/2)/((1/(tau^2))-((pi*w)^2)/4));
t1 = 0:T/1000:T-T/1000; t2 = T:T/1000:2*T-T/1000;
% Vm(0<t<T)=Vm1
Vm1 = Er + c*exp(-t1./tau) + ...
    (I0/(A*Cm))*((1/tau)/((1/(tau^2))-((pi*w)^2)/4))*(sin(pi*w*t1./2)-(pi*w*tau/2)*cos(pi*w*t1./2));
% Vm(t>T)=Vm2
Vm2 = Er + Vm1(end)*exp(-t2./tau);
%final results and graph:
t = [t1,t2]; Vm = [Vm1,Vm2];
figure;
plot(t,Vm); title('Vm progression with Time');
xlabel('Time[msec]'); ylabel('Vm[mV]');
end