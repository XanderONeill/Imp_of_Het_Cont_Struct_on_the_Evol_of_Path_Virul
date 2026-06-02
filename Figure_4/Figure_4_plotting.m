%%% Figure #4 Plotting ESS vs k for different fixed G_T.

clearvars
close all

%% Measures

x = [0.2,0.5,1,2,5,1000000];
y_means = zeros(3, 6);
y_sd = zeros(3, 6);

%% Load data for G_T = 1, find measures

load('Gt1k02.mat');     y_means(1, 1)  = mean(mean_alpha_vec);
                        y_sd(1, 1)     = std(mean_alpha_vec);
load('Gt1k05.mat');     y_means(1, 2)  = mean(mean_alpha_vec);
                        y_sd(1, 2)     = std(mean_alpha_vec);
load('Gt1k1.mat');      y_means(1, 3)  = mean(mean_alpha_vec);
                        y_sd(1, 3)     = std(mean_alpha_vec);
load('Gt1k2.mat');      y_means(1, 4)  = mean(mean_alpha_vec);
                        y_sd(1, 4)     = std(mean_alpha_vec);
load('Gt1k5.mat');      y_means(1, 5)  = mean(mean_alpha_vec);
                        y_sd(1, 5)     = std(mean_alpha_vec);
load('Gt1.mat');        y_means(1, 6)  = mean(mean_alpha_vec);
                       y_sd(1, 6)     = std(mean_alpha_vec);                        

n1 = 1;
p1 = polyfit(x(1:6), y_means(1,:), n1);

%% Load data for G_T = 0, find measures

load('Gt0k02.mat');     y_means(3, 1)  = mean(mean_alpha_vec);
                        y_sd(3, 1)     = std(mean_alpha_vec);
load('Gt0k05.mat');     y_means(3, 2)  = mean(mean_alpha_vec);
                        y_sd(3, 2)     = std(mean_alpha_vec);
load('Gt0k1.mat');      y_means(3, 3)  = mean(mean_alpha_vec);
                        y_sd(3, 3)     = std(mean_alpha_vec);
load('Gt0k2.mat');      y_means(3, 4)  = mean(mean_alpha_vec);
                        y_sd(3, 4)     = std(mean_alpha_vec);
load('Gt0k5.mat');      y_means(3, 5)  = mean(mean_alpha_vec);
                        y_sd(3, 5)     = std(mean_alpha_vec);
load('Gt0.mat');        y_means(3, 6)  = mean(mean_alpha_vec);
                        y_sd(3, 6)     = std(mean_alpha_vec);

% n3 = 1;
% p3 = polyfit(x([1,6]), [y_means(3,1),M3], n3);

%% Figure 1: ESS vs G_T

x = [0.2, 0.5, 1, 2, 5, 100];
y1= y_means - 2*y_sd;
y2= y_means + 2*y_sd;

figure
hold on
% plot(x, y1, 'Color', [0 0.4470 0.7410])
% plot(x, y2, 'Color', [0 0.4470 0.7410])
patch([x fliplr(x)], [y1(1,:) fliplr(y2(1,:))], [0.3010 0.7450 0.9330],'FaceAlpha',.3,'EdgeAlpha',.3)
plot(x, y_means(1,:), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2)
% errorbar(x(1:6), y_means(1,:), 2*y_sd(1,:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
% hold on
% plot(0.2:0.1:10, polyval(p1, 0.2:0.1:10), 'Color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 1.5)

% plot(x, y1, 'Color', [0 0.4470 0.7410])
% plot(x, y2, 'Color', [0 0.4470 0.7410])
patch([x fliplr(x)], [y1(3,:) fliplr(y2(3,:))], [0.6350 0.0780 0.1840],'FaceAlpha',.3,'EdgeAlpha',.3)
plot(x, y_means(3,:), '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
% errorbar(x(1:6), y_means(3,1:6), 2*y_sd(3,1:6), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
% plot(0.2:0.1:10, polyval(p3, 0.2:0.1:10), 'Color', [0.9290 0.6940 0.1250], 'LineStyle', ':', 'LineWidth', 1.5)

xlim([0 10])
ylim([0.49 0.9])
xlabel('Heterogeneity (shape parameter, k)')
ylabel('Evolved level of virulence')

ax = gca;
% set(ax,'xscale','log')
ax.FontSize = 11;
set(gca,'box','off')
%set(gca, 'position', [0.11 0.15 0.36 0.45])
leg1 = legend('G_T = 1', '', 'G_T = 0', 'Location', 'SouthEast');
set(leg1,'Box','off')

grid on
% grid minor
% ax = gca;
% ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

