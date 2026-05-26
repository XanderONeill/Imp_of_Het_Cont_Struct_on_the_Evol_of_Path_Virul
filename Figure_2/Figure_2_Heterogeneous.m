%%% Figure #2 Plotting ESS vs G_T for the heterogeneous (k=0.2) network

clearvars
%close all

% Run this file AFTER the homogeneous run file

%% Measures

x = 0:0.1:1;
y_means = zeros(1, 11);
y_sd = zeros(1, 11);

%% Load data, find measures

load('Gt0k02.mat');     y_means(1)  = mean(mean_alpha_vec(4000:end));
                        y_sd(1)     = std(mean_alpha_vec(end-4000:end));
load('Gt01k02.mat');    y_means(2)  = mean(mean_alpha_vec(4000:end));
                        y_sd(2)     = std(mean_alpha_vec(end-4000:end));
load('Gt02k02.mat');    y_means(3)  = mean(mean_alpha_vec(4000:end));
                        y_sd(3)     = std(mean_alpha_vec(end-4000:end));
load('Gt03k02.mat');    y_means(4)  = mean(mean_alpha_vec(4000:end));
                        y_sd(4)     = std(mean_alpha_vec(end-4000:end));
load('Gt04k02.mat');    y_means(5)  = mean(mean_alpha_vec(4000:end));
                        y_sd(5)     = std(mean_alpha_vec(end-4000:end));
load('Gt05k02.mat');    y_means(6)  = mean(mean_alpha_vec(4000:end));
                        y_sd(6)     = std(mean_alpha_vec(end-4000:end));
load('Gt06k02.mat');    y_means(7)  = mean(mean_alpha_vec(4000:end));
                        y_sd(7)     = std(mean_alpha_vec(end-4000:end));
load('Gt07k02.mat');    y_means(8)  = mean(mean_alpha_vec(4000:end));
                        y_sd(8)     = std(mean_alpha_vec(end-4000:end));
load('Gt08k02.mat');    y_means(9)  = mean(mean_alpha_vec(4000:end));
                        y_sd(9)     = std(mean_alpha_vec(end-4000:end));
load('Gt09k02.mat');    y_means(10) = mean(mean_alpha_vec(4000:end));
                        y_sd(10)    = std(mean_alpha_vec(end-4000:end));
load('Gt1k02.mat');     y_means(11) = mean(mean_alpha_vec(4000:end));
                        y_sd(11)    = std(mean_alpha_vec(end-4000:end));

n = 2;
p = polyfit(x,y_means,n);
p2 = polyfit(x, y_means, n+2);

%% Figure 1: ESS vs G_T

x=0:0.1:1;
y1= y_means - 2*y_sd;
y2= y_means + 2*y_sd;

% plot(x, y1, 'Color', [0 0.4470 0.7410])
% plot(x, y2, 'Color', [0 0.4470 0.7410])
patch([x fliplr(x)], [y1 fliplr(y2)], [0.9290 0.6940 0.1250],'FaceAlpha',.3,'EdgeAlpha',.3)
plot(x, y_means, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
% errorbar(0:0.1:1, y_means, 2*y_sd, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
% hold on
% plot(x, polyval(p, x), 'Color', [0.9290 0.6940 0.1250], 'LineStyle', ':', 'LineWidth', 1.5)
xlim([0 1])
ylim([0.5 1])
yticks([0.5 0.6 0.7 0.8 0.9 1])
xlabel('Proportion of global transmission, G_T')
ylabel('Approximate ESS virulence')

ax = gca;
ax.FontSize = 11;
ylabel('Evolved level of virulence')
set(gca,'box','off')
%set(gca, 'position', [0.11 0.15 0.36 0.45])
% leg1 = legend('k=\infty', '', 'k=0.2', '', 'Location', 'SouthEast');
leg1 = legend('Heterogeneity in Susceptible Contacts', '', 'Heterogeneity in Infected Contacts', '', 'Location', 'SouthEast');
% leg1 = legend('Negative covariance', '', 'Positive covariance', '', 'Location', 'SouthEast');
set(leg1,'Box','off')

grid on
% grid minor
% ax = gca;
% ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

