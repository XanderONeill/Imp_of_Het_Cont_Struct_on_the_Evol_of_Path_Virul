%%% Running the regular network evolution model; 

clearvars
close all
rng("shuffle")


%% Model Set-Up

% GRID SIZE
m = 100; n = 100;     % small worlds grid size is 100 by 100 - this may take some considerable time to run

% ECOLOGICAL AND EPIDEMIOLOGICAL PARAMETERS
b = 10;             % reproduction rate, in Boots & Sasaki this is denoted by r
d = 1;              % natural death rate
alpha_0 = 0.75;      % initial alpha
beta_0 = trade_off_function(alpha_0);
G_T = 0;              % proportion of global transmission
params = [b, d, G_T];

% TIME TO RUN TO
t_end = 16000;
t = 0; t_int = 0;                      % interval start time
t_save = 250;

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

%% Obtaining initial adjacency matrices and digraphs

A_R = RegularAdj(m, n);
G_R = digraph(A_R);
G = G_R;

%% Adding labels to the digraphs

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

degrees_over_time = zeros(length(rpS), t_end+1);
degrees_over_time(:,1) = outdegree(G);

% 1. Mean virulence
mean_alpha_vec = zeros(1, t_end +1);  % keep track of the mean virulence 
                                        % at a given time point
mean_alpha_vec(1) = alpha_0;    % at t = 0, the mean virulence will be 
                                % equal to the initial virulence

%% Stochastic Simulations

tic
while t < t_end
    [G, dt] = Gillespie_Graph(G, G_R, [b, d, G_T]);
    t = t + dt;

    if floor(t) > t_int       % every 5 units of time, we store the mean virulence across the network.
        t_int = floor(t)

        mean_alpha_vec(t_int+1) = mean(G.Nodes.Virulence(G.Nodes.Status=="I"));
        S_over_time(:,t_int+1) = G.Nodes.Status == "S";
        I_over_time(:,t_int+1) = G.Nodes.Status == "I";
        Alpha_over_time(:,t_int+1) = G.Nodes.Virulence;

        degrees_over_time(:, t_int+1) = outdegree(G); 
    end
    if floor(t) >= t_save
        t_save = t_save + 250;
        save('Gt0_InfBirth.mat');
    end
end
toc

save('Gt0_InfBirth')
