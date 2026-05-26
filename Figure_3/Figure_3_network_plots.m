%%% Figure 3: Producing the network visualisation

% Requires the corresponding files which need to be obtained from folder
% Figure_2 in the Github repository.

%% Load in respective file
% Choose one of 'Gt0k02' (heterogeneous, local), 
%               'Gt1k02' (heterogeneous, global),
%               'Gt0' (homogeneous, local),
%               'Gt1' (homogeneous, global) 

% load('Gt0k02');
load('Gt1k02');
% load('Gt0');
% load('Gt1');

O_over_time = ones(size(S_over_time)) - S_over_time - I_over_time;

v = 1:m;
y = repelem(v,m);
x = repmat(1:m, 1, m);

% make all susceptible nodes blue
S_colour = repmat(S_over_time(:,8001),1,3).*[0, 0.4470, 0.7410];
% make all infected nodes red
I_colour = repmat(I_over_time(:,8001),1,3).*[0.6350, 0.0780, 0.1840];
% make all empty nodes white
O_colour = repmat(O_over_time(:,8001),1,3).*[1, 1, 1];
% turn this into node colour data by combining all three
NodeCData = S_colour + I_colour + O_colour;

%% Produce desired figure
fig1 = figure();
p = plot(G, 'XData',x,'YData',y);
p.LineStyle = 'none';
p.MarkerSize = 2;
p.NodeColor = NodeCData;
p.LineWidth = 1.5;
xticklabels('')
yticklabels('')
xticks('');
yticks('');

