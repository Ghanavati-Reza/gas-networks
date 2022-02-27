%% Define Time
dt = 100; % s
tstart = 0;
tend = 20*3600; % s
time = [tstart:dt:tend]';

%% example 1 Definition

% [event time;node index;percentage of change]
t_e = [7*3600 4 0.3; 14*3600 4 0.3];

% number of segments for each pipe
N = [5; 5; 5; 5; 5; 5; 5; 5; 5; 5]; 

% p_s = 12e6 for supply and 
% m_d (mass flow rate) = 400 kg/s for each delivery 

% n = [node: node_index node_type node_spec_value] 
% node_type: 1: supply(pressure spec), 2: junction, 3: delivery(mass flow rate spec)
n = [1 1 12e6; 2 2 0; 3 3 400; 4 3 400; 5 2 0; 6 2 0; 7 2 0; 8 3 400; 9 2 0; 10 2 0; 11 3 400]; 

% p = [pipe: inlet_node_index outlet_node_index] 
p = [1 2; 2 3; 2 4; 2 5; 5 6; 6 7; 7 8; 5 9; 9 10; 10 11]; 

%length of each pipe
L = [10*10^3; 10*10^3; 10*10^3; 10*10^3; 10*10^3; 10*10^3; 10*10^3; 10*10^3; 10*10^3; 10*10^3]; 

%segment length for each pipe
dx = L./N; 

sg = 0.6; % specific gravity of the gas (natural gas)
T = 300; % K
F = 0.0108;
d = 1.0668; % m

%% call ODE solver
[t,x_est,junc] = ode_solver(tstart,tend,t_e,N,p,n,dx,sg,d,F,T);

%% plot the results
plot_results(n,p,N,junc,t,x_est)

