%%% Producing data for the Homogeneous network simulations

clearvars
close all
rng("shuffle")

% Can vary alpha_0 depending on the level of global transmission.
%   For fully homogeneous and 
%       G_T = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
%   use
%   alpha_0 = [0.73, 0.81, 0.9, 0.92, 0.94, 0.96, 0.95, 0.93, 0.92, 0.89, 0.83]

% We do this as it speeds up the simulations if you start nearer to the ESS.

% 1. Set G_T to be what you desire
% 2. Change alpha_0 to be close to the ESS (see above)
% 3. Change save file (final command) to match with correct G_T
% 3. Run - it will take a long time

%% Model Set-Up

% GRID SIZE
m = 100; n = 100;     % run time will grow significantly with an increase in node numbers
                      % would not recommend using anything more than 
                      % 100 x 100

% ECOLOGICAL AND EPIDEMIOLOGICAL PARAMETERS
b = 10;             % reproduction rate, in Boots & Sasaki this is denoted by r
d = 1;              % natural death rate
alpha_0 = 0.92;      % initial alpha
beta_0 = trade_off_function(alpha_0);
G_T = 0.8;              % proportion of global transmission
params = [b, d, G_T];

% TIME TO RUN TO
t_end = 300;
t = 0;
t_int = 0;                      % interval start time

%% Initialisation and setting node status

rpS = zeros(1, m*n); rpI = zeros(1, m*n);

Rand_vec = randperm(m*n, round(0.8*m*n));
rpS(1, Rand_vec(1:round(0.4*m*n))) = 1;         rpI(1, Rand_vec(round(0.4*m*n)+1:round(0.8*m*n))) = 1;
rpO = ones(1, m*n) - rpS - rpI;

node_names = strings(m*n, 1);
for i = 1:n*m
    if rpS(i) == 1
        node_names(i) = "S";
    elseif rpI(i) == 1
        node_names(i) = "I";
    else
        node_names(i) = "O";
    end
end

%% Obtaining initial adjacency matrices

A = RegularAdj(m, n);

%% Creating Birth Graphs and Transmission Graphs

G = digraph(A);

G.Nodes.Status = node_names;

G.Nodes.Virulence(G.Nodes.Status =="I") = alpha_0;
G.Nodes.Transmission(G.Nodes.Status=="I") = beta_0;

%% Data to keep track of

S_over_time = zeros(length(rpS), t_end+1);
S_over_time(:,1) = G.Nodes.Status == "S";
I_over_time = zeros(length(rpI), t_end+1);
I_over_time(:,1) = G.Nodes.Status == "I";
Alpha_over_time = zeros(length(rpI), t_end+1);
Alpha_over_time(:,1) = G.Nodes.Virulence;

% 1. Mean virulence
mean_alpha_vec = zeros(1, t_end +1);  % keep track of the mean virulence 
                                        % at a given time point
mean_alpha_vec(1) = alpha_0;    % at t = 0, the mean virulence will be 
                                % equal to the initial virulence
% 2. Shelf Shading ability
Self_shade = zeros(1, t_end+1);   % keep track of the total number of 
                                    % local susceptibles around all
                                    % infected nodes
Self_shade(1) = (G.Nodes.Status=="I")'*((G.Nodes.Status=="S")'*G.adjacency')';


%% Stochastic Simulations


while t < t_end
    [G, dt] = Gillespie_Graph(G, [b, d, G_T]);
    t = t + dt;

    if floor(t) > t_int       % every 5 units of time, we store the mean virulence across the network.
        t_int = floor(t)

        mean_alpha_vec(t_int+1) = mean(G.Nodes.Virulence(G.Nodes.Status=="I"));
        S_over_time(:,t_int+1) = G.Nodes.Status == "S";
        I_over_time(:,t_int+1) = G.Nodes.Status == "I";
        Alpha_over_time(:,t_int+1) = G.Nodes.Virulence;
    end
end

save('Gt08')
