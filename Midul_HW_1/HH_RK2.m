function [solution,nrk,mrk,hrk] = HH_RK2(dt, Vm, n, m, h, I)
% Define the needed constants
Gna = 1.2;  % mS/mm^2
Gk = 0.36;  % mS/mm^2
Gl = 0.003; % mS/mm^2
Ena = 50;   % mV
Ek = -77;   % mV
El = -54.4; % mV
C = 0.01;   % uF/nm^2
V = Vm;
nrk = zeros(length(Vm),1);
mrk = zeros(length(Vm),1);
hrk = zeros(length(Vm),1);
nrk(1)=n;mrk(1)=m;hrk(1)=h;
for i = 2:length(Vm)

    % Parameters needed for RK step
    tau = [calc_gate('tau_n',V(i-1)),calc_gate('tau_m',V(i-1)),calc_gate('tau_h',V(i-1))];
    inf = [calc_gate('n_inf',V(i-1)),calc_gate('m_inf',V(i-1)),calc_gate('h_inf',V(i-1))];
    
    % Derivatives
    dgate_dt = (inf-[nrk(i-1),mrk(i-1),hrk(i-1)])./(tau);
    dV_dt = Voltage_derivative(C, Gk, Gna, Gl, nrk(i-1), mrk(i-1), hrk(i-1), I(i-1),...
        V(i-1), Ek, Ena, El);

    % K1 = h*f(t,V)
    K1 = dt.*[dgate_dt, dV_dt];
    

    % K2 = h*f(t+0.5h,V+0.5*K1)
    tau = [calc_gate('tau_n',V(i-1)+0.5*K1(1)),calc_gate(...
        'tau_m',V(i-1)+0.5*K1(2)),calc_gate('tau_h',V(i-1)+0.5*K1(3))];
    inf = [calc_gate('n_inf',V(i-1)+0.5*K1(1)),calc_gate(...
        'm_inf',V(i-1)+0.5*K1(2)),calc_gate('h_inf',V(i-1)+0.5*K1(3))];
    
    % Derivatives
    dgate_dt = dt*(inf-[nrk(i-1),mrk(i-1),hrk(i-1)])./(tau);
    dV_dt = dt.*Voltage_derivative(C,Gk,Gna,Gl,nrk(i-1),mrk(i-1),hrk(i-1),...
        I(i-1),V(i-1)+0.5*K1(4),Ek,Ena,El);

    % Gates and Voltage update after RK step
    nrk(i) = nrk(i-1) + dgate_dt(1);
    mrk(i) = mrk(i-1) + dgate_dt(2);
    hrk(i) = hrk(i-1) + dgate_dt(3);
    V(i) = V(i-1) + dV_dt;

end

solution = V;
end
%% Local Functions
% Volatage derivative
function [dV_dt] = Voltage_derivative(C, Gk, Gna, Gl, n, m, h, I, Vm, Ek, Ena, El)
    dV_dt = (I - Gk*(n^4)*(Vm-Ek) - Gna*(m^3)*h*(Vm - Ena) - Gl*(Vm-El))/(C);
end

% inf and tau calculations
function [gate]=calc_gate(type,v)
    if strcmp(type,'n_inf')
        gate = calc_a_n(v)./(calc_a_n(v) + calc_b_n(v));
    elseif strcmp(type,'m_inf')
        gate = calc_a_m(v)./(calc_a_m(v) + calc_b_m(v));
    elseif strcmp(type,'h_inf')
        gate = calc_a_h(v)./(calc_a_h(v) + calc_b_h(v));
    elseif strcmp(type,'tau_n')
        gate = 1./(calc_a_n(v) + calc_b_n(v));
    elseif strcmp(type,'tau_m')
        gate = 1./(calc_a_m(v) + calc_b_m(v));
    elseif strcmp(type,'tau_h')
        gate = 1./(calc_a_h(v) + calc_b_h(v));
    end
end

% a and b calculations for inf and tau
function a_n = calc_a_n(v)
    a_n = (0.01*(v+55))./(1-exp(-0.1*(v+55)));
end

function a_m = calc_a_m(v)
    a_m = (0.1*(v+40))./(1-exp(-0.1*(v+40)));
end

function a_h = calc_a_h(v)
    a_h = 0.07*exp(-0.05*(v+65));
end

function b_n = calc_b_n(v)
    b_n = 0.125*exp(-0.0125*(v+65));
end

function b_m = calc_b_m(v)
    b_m = 4*exp(-0.0556*(v+65));
end

function b_h = calc_b_h(v)
    b_h = 1./(1 + exp(-0.1*(v+35)));
end
