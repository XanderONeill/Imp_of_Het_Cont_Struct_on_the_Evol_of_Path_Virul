%%% Figure #2 Plotting ESS vs G_T for the homogeneous network

clearvars
close all

%%% Run this file before the heterogeneous network file

%% Measures

x = 0:0.1:1;
y_means = zeros(1, 11);
y_sd = zeros(1, 11);

%% Load data, find measures

load('Gt0.mat');    y_means(1)  = mean(mean_alpha_vec(4000:end));
                    y_sd(1)     = std(mean_alpha_vec(end-4000:end));
load('Gt01.mat');   y_means(2)  = mean(mean_alpha_vec(4000:end));
                    y_sd(2)     = std(mean_alpha_vec(end-4000:end));
load('Gt02.mat');   y_means(3)  = mean(mean_alpha_vec(4000:end));
                    y_sd(3)     = std(mean_alpha_vec(end-4000:end));
load('Gt03.mat');   y_means(4)  = mean(mean_alpha_vec(4000:end));
                    y_sd(4)     = std(mean_alpha_vec(end-4000:end));
load('Gt04.mat');   y_means(5)  = mean(mean_alpha_vec(4000:end));
                    y_sd(5)     = std(mean_alpha_vec(end-4000:end));
load('Gt05.mat');   y_means(6)  = mean(mean_alpha_vec(4000:end));
                    y_sd(6)     = std(mean_alpha_vec(end-4000:end));
load('Gt06.mat');   y_means(7)  = mean(mean_alpha_vec(4000:end));
                    y_sd(7)     = std(mean_alpha_vec(end-4000:end));
load('Gt07.mat');   y_means(8)  = mean(mean_alpha_vec(4000:end));
                    y_sd(8)     = std(mean_alpha_vec(end-4000:end));
load('Gt08.mat');   y_means(9)  = mean(mean_alpha_vec(4000:end));
                    y_sd(9)     = std(mean_alpha_vec(end-4000:end));
load('Gt09.mat');   y_means(10) = mean(mean_alpha_vec(4000:end));
                    y_sd(10)    = std(mean_alpha_vec(end-4000:end));
load('Gt1.mat');    y_means(11) = mean(mean_alpha_vec(4000:end));
                    y_sd(11)    = std(mean_alpha_vec(end-4000:end));

n = 2;
p = polyfit(x,y_means,n);
p1 = polyfit(x, y_means, n+1);

%% Figure 1: ESS vs G_T

x=0:0.1:1;
y1= y_means - 2*y_sd;
y2= y_means + 2*y_sd;

fig = figure;
hold on
% plot(x, y1, 'Color', [0 0.4470 0.7410])
% plot(x, y2, 'Color', [0 0.4470 0.7410])
patch([x fliplr(x)], [y1 fliplr(y2)], [0 0.4470 0.7410],'FaceAlpha',.3,'EdgeAlpha',.3)
plot(x, y_means, '-.','Color', [0 0.4470 0.7410], 'LineWidth', 2)
% errorbar(0:0.1:1, y_means, 2*y_sd, 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
% hold on
% plot(x, polyval(p, x), 'Color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 1.5)
xlim([0 1])
xlabel('Proportion of global transmission, G_T')
ylim([0.5 1])



