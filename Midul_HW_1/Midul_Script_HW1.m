%% Midul HW1 Script
% pay attention to which section represent each Hw section
% make sure all mentioned functions are in the same file as script
clear all; clc; close all;
%% Q2b
% show graph for Vm computed in Q2a, which constants defined here:
% set parameters
I0 = 10^-9;         % A
tau = 5;            % msec
T = 800;            % msec
Cm = 1;             % nF/mm^2
R = 5*10^-6;        % m
Er = -65;           % mV
% make sure 'graphing_Q1' is in file
[~] = graphing_Q2(Er,I0,T,R,Cm,tau);
%% Q3
% define parameters for the problem
% other params as E,g are define within the function
dt = 0.01;                                   % msec
tmesh = 0:dt:1000-dt;                        % msec
Vm = [-60, zeros(length(tmesh)-1,1)'];       % mV
h = 0; n = 0; m = 0;

%% Q3a
% I = 0 nA
% make sure 'HH_RK2' is in file
I = zeros(length(tmesh),1); % nA
[sol,n_new,m_new,h_new] = HH_RK2(dt, Vm, n, m, h, I);
figure;
plot(tmesh,sol);hold on; title('IV curve for I=0');
plot(tmesh,I); xlabel('Time[msec]'); ylabel('Voltage or Current[mV or nA]');
legend('Voltage','Current'); hold off;
figure;
plot(tmesh,n_new);hold on; title('n,m,h for I=0');
plot(tmesh,m_new); plot(tmesh,h_new);
xlabel('Time[msec]'); ylabel('Probability');
% xlim([0 50]);
legend('n gates','m gates','h gates'); hold off;
%% Q3b
% I = 1 nA (250msec < time < 750msec)
I = [zeros(25000,1);ones(50000,1);zeros(25000,1)]; % nA
[sol,~] = HH_RK2(dt, Vm, n, m, h, I);
figure;
plot(tmesh,sol);hold on; title('IV curve for I=1nA for 500 msec');
plot(tmesh,I); xlabel('Time[msec]'); ylabel('Voltage or Current[mV or nA]');
legend('Voltage','Current'); hold off;
%% Q3c
% I = 1 nA (250msec < time < 750msec)
% m,n,h in accordance to estimated values of newron in Resting State
% Stabilized values from Q3a
n = n_new(end); m = m_new(end); h = h_new(end);
[sol,~] = HH_RK2(dt, Vm, n, m, h, I);
figure;
plot(tmesh,sol);hold on; title('IV curve for I=1nA for 500 msec, Resting State gates probabilities');
plot(tmesh,I); xlabel('Time[msec]'); ylabel('Voltage or Current[mV or nA]');
legend('Voltage','Current'); hold off;
%% Q3d
% I = 0.05 + 0.1*sin(20*pi*t)
% m,n,h in accordance to estimated values of newron in Resting State
% Stabilized values from Q3a
sin_vec = 0.05 + 0.1*sin(20*pi*tmesh(25000:74999)/1000); sin_vec=sin_vec';
I = [zeros(25000,1);sin_vec;zeros(25000,1)]; % nA
[sol,~] = HH_RK2(dt, Vm, n, m, h, I);
figure;
subplot(2,1,1);
plot(tmesh,sol); title('V for I = 0.05+0.1sin(20*pi*t) nA');
xlabel('Time[msec]'); ylabel('Voltage[mV]');
subplot(2,1,2);
plot(tmesh,I);  title('I = 0.05+0.1sin(20*pi*t) nA');
xlabel('Time[msec]'); ylabel('Current[nA]');




