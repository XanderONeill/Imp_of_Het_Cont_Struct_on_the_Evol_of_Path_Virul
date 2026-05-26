%%% Producing data for Figure 4
% where we fix the globality/locality of transmission and vary the
% heterogeneity levels

clearvars
close all
rng("shuffle")

% Requires functions: (i) RegularAdj, (ii) IrregularAdj, 
% (iii) Gillespie_Graph, and (iv) trade_off_function

% 1. Set G_T to be what you desire
% 2. Set k to be what you desire
% 3. Change alpha_0 to be close to the ESS (see Figure 4 if using similar values)
% 4. Change save file (final command) to match with correct G_T
% 5. Run - it will take a long time

% Note if you have already produced Figure 2, you will likely have some of
% the save data files already.

%% Model Set-Up

% GRID SIZE
m = 100; n = 100;     % small worlds grid size is 100 by 100 - this may take some considerable time to run

% GAMMA DISTRIBUTION
k = 0.2;              % k -> infinity means fully homogeneous
mu = 4;             % average number of neighbours 

% ECOLOGICAL AND EPIDEMIOLOGICAL PARAMETERS
b = 10;             % reproduction rate, in Boots & Sasaki this is denoted by r
d = 1;              % natural death rate
alpha_0 = 0.56;      % initial alpha
beta_0 = trade_off_function(alpha_0);
G_T = 0;              % proportion of global transmission
params = [b, d, G_T];

% TIME TO RUN TO
t_end = 1000;
t = 0;
t_int = 0;                      % interval start time

%% Initialisation and setting node status

rpS = zeros(1, m*n); rpI = zeros(1, m*n);

Rand_vec = randperm(m*n, round(0.6*m*n));
rpS(1, Rand_vec(1:round(0.2*m*n))) = 1;         rpI(1, Rand_vec(round(0.2*m*n)+1:round(0.6*m*n))) = 1;
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

A_R = RegularAdj(m, n);
A_I = IrregularAdj(m, n, k, mu);
%save('A_I_k02')

%% Creating Birth Graphs and Transmission Graphs

G_R = digraph(A_R);
G = digraph(A_I);

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
    [G, dt] = Gillespie_Graph(G, G_R, [b, d, G_T]);
    t = t + dt;

    if floor(t) > t_int       % every 5 units of time, we store the mean virulence across the network.
        t_int = floor(t)

        mean_alpha_vec(t_int+1) = mean(G.Nodes.Virulence(G.Nodes.Status=="I"));
        S_over_time(:,t_int+1) = G.Nodes.Status == "S";
        I_over_time(:,t_int+1) = G.Nodes.Status == "I";
        Alpha_over_time(:,t_int+1) = G.Nodes.Virulence;
        
        Self_shade(t_int+1) = (G.Nodes.Status=="I")'*((G.Nodes.Status=="S")'*G.adjacency')';
    end
end

save('Gt0k02')
