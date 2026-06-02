%%% Script to produce Figure 5 and Figure S2 in the paper

% Requires the running of scripts:
% Hetero_birth
% Infect_birth
% Global_reprod
% to gain additional data files
% Files 'Gt0.mat' and 'Gt0k02' can be copied over from previous simulations

means = zeros(8,1);
stds = zeros(8,1);

%% Original Model Results
load('Gt0.mat')
means(1,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(1,1) = std(mean_alpha_vec(t_int-4000:t_int));

load('Gt0k02.mat')
means(2,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(2,1) = std(mean_alpha_vec(t_int-4000:t_int));

%% Heterogeneity at birth results
load('Gt0.mat')
means(3,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(3,1) = std(mean_alpha_vec(t_int-4000:t_int));

load('Gt0k02_birth.mat')
means(4,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(4,1) = std(mean_alpha_vec(t_int-4000:t_int));

%% Infected Give Birth Results
load('Gt0_InfBirth.mat')
means(5,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(5,1) = std(mean_alpha_vec(t_int-4000:t_int));

load('Gt0k02_InfBirth_ext.mat')
means(6,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(6,1) = std(mean_alpha_vec(t_int-4000:t_int));

%% Global Reproduction results
load('Gt0_GR.mat')
means(7,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(7,1) = std(mean_alpha_vec(t_int-4000:t_int));

load('Gt0k02_GR_ext.mat')
means(8,1) = mean(mean_alpha_vec(t_int-4000:t_int));
stds(8,1) = std(mean_alpha_vec(t_int-4000:t_int));

%% Plotting the results
close all
x = 1:4;
y1 = means(2*x-1,:);
y2 = means(2*x,:);
err1 = stds(2*x-1,:);
err2 = stds(2*x,:);

%% Figure S1: Heterogeneity at birth vs original
close all
fig1 = figure;
e1 = errorbar(x(1:2), y1(1:2), 2*err1(1:2), 'o', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on
e2 = errorbar(x(1:2), y2(1:2), 2*err2(1:2), 'o', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250]);

e1.MarkerFaceColor = [0.2 0.6470 0.9410];
e1.MarkerSize = 5;

e2.MarkerFaceColor = [0.95 0.75 0.15];
e2.MarkerSize = 5;

xlim([0.5 2.5])
xticks([1, 2])
xticklabels({'Original', 'Heterogeneity at birth'})
xtickangle(45)

ylabel('Evolved level of virulence')
set(gca, 'Position', [0.15 0.25 0.5 0.65])

leg1 = legend('Regular, G_T = 0', 'Irregular, G_T = 0', 'Position', [0.75 0.7 0.05 0.05]);
legend boxoff;

%% Figure 5: Global reproduction and infecteds give birth

fig2 = figure;
e1 = errorbar(x(1:3), y1([1, 3:4]), 2*err1([1, 3:4]), 'o', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
hold on
e2 = errorbar(x(1:3), y2([1, 3:4]), 2*err2([1, 3:4]), 'o', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250]);

e1.MarkerFaceColor = [0.2 0.6470 0.9410];
e1.MarkerSize = 5;

e2.MarkerFaceColor = [0.95 0.75 0.15];
e2.MarkerSize = 5;

xlim([0.75 3.25])
xticks([1, 2, 3])
xticklabels({'Original', 'All individuals give birth', 'Global reproduction'})
xtickangle(45)
set(gca, 'Position', [0.15 0.3 0.5 0.65])

ylabel('Evolved level of virulence')

leg2 = legend('Regular, G_T = 0', 'Irregular, G_T = 0', 'Position', [0.75 0.75 0.05 0.05]);
legend boxoff;