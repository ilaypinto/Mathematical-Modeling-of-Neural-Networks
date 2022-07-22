function [solution,Psrk,Zrk] = LIF_RK2_2N(dt, Vm, Es)
El = -70;
tau_m = 20;
tau_s = 10;
thresh = -54;
Pmax = 0.5;
rmgs = 0.15;
RmIe = 18;
Zrk = zeros(length(0:dt:5000-dt),2);
Psrk = zeros(length(0:dt:5000-dt),2);
V = zeros(length(0:dt:5000-dt),2);
V(1,1) = Vm(1);
V(1,2) = Vm(2);
    for i = 2:size(V,1)
        dZ_dt = zeros(2,1); dPs_dt = zeros(2,1); dV_dt = zeros(2,1);

        % First neuron
        dZ_dt(1) = Z_derivative(Zrk(i-1,1),tau_s);
        dPs_dt(1) = Ps_derivative(Pmax,Zrk(i-1,1),Psrk(i-1,1),tau_s);
        dV_dt(1) = Voltage_derivative(V(i-1,1),tau_m,El,RmIe,rmgs,Es,Psrk(i-1,1));

        % Second neuron
        dZ_dt(2) = Z_derivative(Zrk(i-1,2),tau_s);
        dPs_dt(2) = Ps_derivative(Pmax,Zrk(i-1,2),Psrk(i-1,2),tau_s);
        dV_dt(2) = Voltage_derivative(V(i-1,2),tau_m,El,RmIe,rmgs,Es,Psrk(i-1,2));

        % K1 = h*f(t,V)
        K1_1 = dt.*[dZ_dt(1),dPs_dt(1),dV_dt(1)];
        K1_2 = dt.*[dZ_dt(2),dPs_dt(2),dV_dt(2)];
        
        % K2 = h*f(t+0.5h,V+0.5*K1)
        dZ_dt(1) = dt.*Z_derivative(Zrk(i-1,1)+0.5*K1_1(1),tau_s);
        dPs_dt(1) = dt.*Ps_derivative(Pmax,Zrk(i-1,1)+0.5*K1_1(1),Psrk(i-1,1)+0.5*K1_1(2),tau_s);
        dV_dt(1) = dt.*Voltage_derivative(V(i-1,1)+0.5*K1_1(3),tau_m,El,RmIe,rmgs,Es,Psrk(i-1,1)+0.5*K1_1(2));

        dZ_dt(2) = dt.*Z_derivative(Zrk(i-1,2)+0.5*K1_2(1),tau_s);
        dPs_dt(2) = dt.*Ps_derivative(Pmax,Zrk(i-1,2)+0.5*K1_2(1),Psrk(i-1,2)+0.5*K1_2(2),tau_s);
        dV_dt(2) = dt.*Voltage_derivative(V(i-1,2)+0.5*K1_2(3),tau_m,El,RmIe,rmgs,Es,Psrk(i-1,2)+0.5*K1_2(2));
    
        % Gates and Voltage update after RK step
        Zrk(i,1) = Zrk(i-1,1) + dZ_dt(1);
        Psrk(i,1) = Psrk(i-1,1) + dPs_dt(1);
        V(i,1) = V(i-1,1) + dV_dt(1);

        Zrk(i,2) = Zrk(i-1,2) + dZ_dt(2);
        Psrk(i,2) = Psrk(i-1,2) + dPs_dt(2);
        V(i,2) = V(i-1,2) + dV_dt(2);

        if V(i,1) > thresh
            V(i,1) = V(1,1);
            Zrk(i,2) = 1;
        end
        if V(i,2) > thresh
            V(i,2) = V(1,1);
            Zrk(i,1) = 1;
        end
    
    end

    solution = V;
end

%% Local Functions
% Volatage derivative
function [dV_dt] = Voltage_derivative(Vm,tau_m,El,RmIe,rmgs,Es,Ps)
    dV_dt = (El - Vm + RmIe - rmgs*Ps*(Vm - Es))/(tau_m);
end

function [dPs_dt] = Ps_derivative(Pmax,Z,Ps,tau_s)
    dPs_dt = (exp(1)*Pmax*Z*(1-Ps) - Ps)/(tau_s);
end

function [dZ_dt] = Z_derivative(Z,tau_s)
    dZ_dt = -1*Z/tau_s;
end