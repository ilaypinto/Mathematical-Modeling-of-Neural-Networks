function [solution,spikes,Psrk,Zrk] = LIF_RK2(dt, Vm, section,synapse)
spikes = []; 
El = -70;
tau_m = 20;
tau_s = 10;
thresh = -54;

if section == 1

% Section 1- No synapse involved

    rmgs = 0;
    Es = 0;
    Ps = 0;
    Pmax = 0;
    Psrk = 0;
    Zrk = 0;
    RmIe = 25;
    V = zeros(length(0:dt:500-dt),1);
    V(1) = Vm;
    
    for i = 2:length(V)
    
        dV_dt = Voltage_derivative(V(i-1),tau_m,El,RmIe,rmgs,Es,Ps);
    
        % K1 = h*f(t,V)
        K1 = dt.*dV_dt;
        
        % K2 = h*f(t+0.5h,V+0.5*K1)
        dV_dt = dt.*Voltage_derivative(V(i-1)+0.5*K1,tau_m,El,RmIe,rmgs,Es,Ps);
    
        % Gates and Voltage update after RK step
        V(i) = V(i-1) + dV_dt;
        if V(i) > thresh
            V(i) = V(1);
            spikes = [spikes, i];
        end
    
    end

    solution = V;
    spikes = spikes/100;

elseif section == 2

% Excitatory synapse, worse model for Ps

    rmgs = 0.5;
    Es = 0;
    RmIe = 15;
    Pmax = 0.5;  % Change to Pmax = 1 before running section 3 on script
    Zrk = zeros(length(0:dt:500-dt),1);
    Psrk = zeros(length(0:dt:500-dt),1);
    V = zeros(length(0:dt:500-dt),1);
    V(1) = Vm;
    for i = 2:length(V)

        dZ_dt = Z_derivative(Zrk(i-1),tau_s);
        dPs_dt = Ps_derivative2(Pmax,Zrk(i-1),Psrk(i-1),tau_s);
        dV_dt = Voltage_derivative(V(i-1),tau_m,El,RmIe,rmgs,Es,Psrk(i-1));
    
        % K1 = h*f(t,V)
        K1 = dt.*[dZ_dt,dPs_dt,dV_dt];
        
        % K2 = h*f(t+0.5h,V+0.5*K1)
        dZ_dt = dt.*Z_derivative(Zrk(i-1)+0.5*K1(1),tau_s);
        dPs_dt = dt.*Ps_derivative2(Pmax,Zrk(i-1)+0.5*K1(1),Psrk(i-1)+0.5*K1(2),tau_s);
        dV_dt = dt.*Voltage_derivative(V(i-1)+0.5*K1(3),tau_m,El,RmIe,rmgs,Es,Psrk(i-1)+0.5*K1(2));
    
        % Gates and Voltage update after RK step
        Zrk(i) = Zrk(i-1) + dZ_dt;
        Psrk(i) = Psrk(i-1) + dPs_dt;
        V(i) = V(i-1) + dV_dt;
        if ismember(i,synapse*100)
             Zrk(i) = 1;
        end
        if V(i) > thresh
            V(i) = V(1);
            spikes = [spikes, i];
        end
    
    end

    solution = V;
    spikes = spikes/100;

elseif section == 3

% Excitatory synapse, better model for Ps

    rmgs = 0.5;
    Es = 0;
    RmIe = 15;
    Pmax = 1;
    Zrk = zeros(length(0:dt:500-dt),1);
    Psrk = zeros(length(0:dt:500-dt),1);
    V = zeros(length(0:dt:500-dt),1);
    V(1) = Vm;
    for i = 2:length(V)

        dZ_dt = Z_derivative(Zrk(i-1),tau_s);
        dPs_dt = Ps_derivative3(Pmax,Zrk(i-1),Psrk(i-1),tau_s);
        dV_dt = Voltage_derivative(V(i-1),tau_m,El,RmIe,rmgs,Es,Psrk(i-1));
    
        % K1 = h*f(t,V)
        K1 = dt.*[dZ_dt,dPs_dt,dV_dt];
        
        % K2 = h*f(t+0.5h,V+0.5*K1)
        dZ_dt = dt.*Z_derivative(Zrk(i-1)+0.5*K1(1),tau_s);
        dPs_dt = dt.*Ps_derivative3(Pmax,Zrk(i-1)+0.5*K1(1),Psrk(i-1)+0.5*K1(2),tau_s);
        dV_dt = dt.*Voltage_derivative(V(i-1)+0.5*K1(3),tau_m,El,RmIe,rmgs,Es,Psrk(i-1)+0.5*K1(2));
    
        % Gates and Voltage update after RK step
        Zrk(i) = Zrk(i-1) + dZ_dt;
        Psrk(i) = Psrk(i-1) + dPs_dt;
        V(i) = V(i-1) + dV_dt;
        if ismember(i,synapse*100)
             Zrk(i) = 1;
        end
        if V(i) > thresh
            V(i) = V(1);
            spikes = [spikes, i];
        end
    
    end

    solution = V;
    spikes = spikes/100;
end
end
%% Local Functions
% Volatage derivative
function [dV_dt] = Voltage_derivative(Vm,tau_m,El,RmIe,rmgs,Es,Ps)
    dV_dt = (El - Vm + RmIe - rmgs*Ps*(Vm - Es))/(tau_m);
end

function [dPs_dt] = Ps_derivative2(Pmax,Z,Ps,tau_s)
    dPs_dt = (exp(1)*Pmax*Z - Ps)/(tau_s);
end

function [dPs_dt] = Ps_derivative3(Pmax,Z,Ps,tau_s)
    dPs_dt = (exp(1)*Pmax*Z*(1-Ps) - Ps)/(tau_s);
end

function [dZ_dt] = Z_derivative(Z,tau_s)
    dZ_dt = -1*Z/tau_s;
end
