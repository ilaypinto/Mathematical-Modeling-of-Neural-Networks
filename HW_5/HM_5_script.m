%% Midul HW 5 Script %%

clear all; clc;close all;
%% Section 1

% LIF model with no Synapses
figure;
[solution,spikes,~,~] = LIF_RK2(0.01,-80, 1,0);
plot(0:0.01:500-0.01,solution); title('LIF model with no Synapses');
xline(spikes); xlabel('Time[msec]');ylabel('Voltage[mV]')


%% Section 2

% LIF model with Excitatory Synapse
figure;
synapse = [50,60,150,190,300,320,400,410,450];
[solution2,spikes,Psrk,~] = LIF_RK2(0.01,-80, 2,synapse);
subplot(2,1,1);
plot(0:0.01:500-0.01,solution2); title('LIF model with Excitatory Synapse');
xline(spikes);xlabel('Time[msec]');ylabel('Voltage[mV]')
subplot(2,1,2);
plot(0:0.01:500-0.01,Psrk);
xlabel('Time[msec]');ylabel('Probability')

%% Section 3

% Comparison of Ps models
% Before running this section, change Pmax in LIF_RK2 in section 2 to
% Pmax = 1

% 1Hz
figure;
synapse = 1;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 2,synapse);
plot(0:0.01:500-0.01,Prk); title('Comparison of Ps models: 1Hz');
hold on;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 3,synapse);
plot(0:0.01:500-0.01,Prk); xlabel('Time[msec]');ylabel('Probability');
legend({'dPs_dt=(exp(1)*Pmax*Z - Ps)/(tau_s)','dPs_dt=(exp(1)*Pmax*Z*(1-Ps) - Ps)/(tau_s)'},'Location','southeast');
hold off;

% 10Hz
figure;
synapse =0.5:100:499.5;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 2,synapse);
plot(0:0.01:500-0.01,Prk); title('Comparison of Ps models: 10Hz');
hold on;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 3,synapse);
plot(0:0.01:500-0.01,Prk);xlabel('Time[msec]');ylabel('Probability');
legend({'dPs_dt=(exp(1)*Pmax*Z - Ps)/(tau_s)','dPs_dt=(exp(1)*Pmax*Z*(1-Ps) - Ps)/(tau_s)'},'Location','southeast')
hold off;

% 50Hz
figure;
synapse =0.5:10:499.5;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 2,synapse);
plot(0:0.01:500-0.01,Prk); title('Comparison of Ps models: 50Hz');
hold on;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 3,synapse);
plot(0:0.01:500-0.01,Prk);xlabel('Time[msec]');ylabel('Probability');
legend({'dPs_dt=(exp(1)*Pmax*Z - Ps)/(tau_s)','dPs_dt=(exp(1)*Pmax*Z*(1-Ps) - Ps)/(tau_s)'},'Location','southeast')
hold off;

% 100Hz
figure;
synapse =0.5:5:499.5;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 2,synapse);
plot(0:0.01:500-0.01,Prk); title('Comparison of Ps models: 100Hz');
hold on;
[~,~,Prk,~] = LIF_RK2(0.01,-80, 3,synapse);
plot(0:0.01:500-0.01,Prk);xlabel('Time[msec]');ylabel('Probability');
legend({'dPs_dt=(exp(1)*Pmax*Z - Ps)/(tau_s)','dPs_dt=(exp(1)*Pmax*Z*(1-Ps) - Ps)/(tau_s)'},'Location','southeast')
hold off;

%% Section 4

% Two connected neurons: Excitatory synapse
[solution,~,~] = LIF_RK2_2N(0.01,[-80,-70], 0);
figure;
plot(0:0.01:5000-0.01,solution(:,1)); title('Two connected neurons: Excitatory synapse');
hold on;
plot(0:0.01:5000-0.01,solution(:,2));xlabel('Time[msec]');ylabel('Voltage[mV]')
legend({'V0 = -80mV','V0 = -70mV'},'Location','southeast')
hold off;

%% Section 5

% Two connected neurons: Inhibitory synapse
[solution,~,~] = LIF_RK2_2N(0.01,[-80,-70], -80);
figure;
plot(0:0.01:5000-0.01,solution(:,1)); title('Two connected neurons: Inhibitory synapse');
hold on;
plot(0:0.01:5000-0.01,solution(:,2));xlabel('Time[msec]');ylabel('Voltage[mV]') 
legend({'V0 = -80mV','V0 = -70mV'},'Location','southeast')
hold off;