clearvars
close all
rng("shuffle")

% Code that can reproduce Figure 1 in the main manuscript
% Use k = 0.2 for heterogeneous case
% Use k = large number for homogeneous case (e.g. k = 10000)

% Note, repeatedly running this script for k = 100000 will reproduce the
% exact same plots as in the manuscript. However, k = 0.2 will not due to
% the randomness of the initialisation.

% Uses the functions RegularAdj(m,n) and IrregularAdj(G_R, k, mu), which
% themselves produce homogeneous and heterogeneous adjacency matrices,
% respectively. 

%% Model Set-Up

% GRID SIZE
m = 100; n = 100;     % small worlds grid size is 100 by 100 - this may take some considerable time to run

% GAMMA DISTRIBUTION
k = 0.2;              % k -> infinity means fully homogeneous
mu = 4;             % average number of neighbours 

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

A_I = IrregularAdj(G_R, k, mu);
G = digraph(A_I);

%% Adding labels to the digraphs

G.Nodes.Status = node_names;

v = 1:m;
y = repelem(v,m);
x = repmat(1:m, 1, m);

%% Figure 1, Network Subfigure

subfig1 = figure();
id = [441 442 443 444 445 541 542 543 544 545 641 642 643 644 645 741 742 743 744 745 841 842 843 844 845];
H = subgraph(G,id);

idx = repelem(1:5, 5);
idy = repmat(1:5, 1, 5);

p1 = plot(H,'XData',idx,'YData',idy, 'NodeLabel', '');

H_NodeCData = zeros(numnodes(H), 3);
for i = 1:numnodes(H)
    if G.Nodes.Status(id(i))=='S'
        H_NodeCData(i,:) = [0 0.4470 0.7410];
    elseif G.Nodes.Status(id(i)) == 'I'
        H_NodeCData(i,:) = [0.6350 0.0780 0.1840];
    else
        H_NodeCData(i,:) = [.8 .8 .8];
    end
end

p1.LineStyle = '--';
p1.MarkerSize = 7;
p1.NodeColor = H_NodeCData;
p1.LineWidth = 1.5;
p1.ArrowSize = 10;
xticklabels('')
yticklabels('')
xticks('');
yticks('');

%% Figure 1, Histogram subfigure

subfigure2 = figure();
h = histogram(outdegree(G), -0.5:1:30.5);
set(gca, 'Position', [0.15 0.15 0.75 0.75])
set(gca, 'FontSize', 14)
ylabel('Number of nodes with specific degree');
xlabel('Individuals outgoing degree');

%% Figure 3 style network plot, pre epi and evo dynamics

% p = plot(G, 'XData',x,'YData',y);% 'NodeLabel', G.Nodes.Status);
% p.LineStyle = 'none';
% p.MarkerSize = 2;
% 
% NodeCData = zeros(m*n, 3);
% for i = 1:m*n
%     if G.Nodes.Status(i)=='S'
%         NodeCData(i,:) = [0 0.4470 0.7410];
%     elseif G.Nodes.Status(i) == 'I'
%         NodeCData(i,:) = [0.6350 0.0780 0.1840];
%     else
%         NodeCData(i,:) = [1 1 1];
%     end
% end
% 
% p.NodeColor = NodeCData;
% p.LineWidth = 1.5;
% xticklabels('')
% yticklabels('')
% xticks('');
% yticks('');