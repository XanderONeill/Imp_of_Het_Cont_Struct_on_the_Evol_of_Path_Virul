clearvars
close all

% Requires the obtaining of data from the homogeneous network simulations.
% Cannot upload our data files to Github as they take up too much file
% space, but make sure to run simulations to obtain the correct data files
% first.

%% Assortativity of any individual
% We ignore empty nodes as no individuals are present

R_over_global_over_time = zeros(11, 4001);  % assortativity of individuals as a whole

%% Assortativity of infected individuals only
II_over_global_over_time =  zeros(11, 4001);              

%% Assortativity of susceptible individuals only
SS_over_global_over_time =  zeros(11, 4001);              

%% Load relevant files

S = zeros(10000, 4001, 11); I = zeros(10000, 4001, 11); O = zeros(10000, 4001, 11);
% Regular networks
load('Gt0');    S(:,:,1) = S_over_time(:,end-4000:end); 
                I(:,:,1) = I_over_time(:,end-4000:end); 
                O(:,:,1) = ones(10000, 4001, 1) - S(:,:,1) - I(:,:,1);
load('Gt01');   S(:,:,2) = S_over_time(:,end-4000:end); 
                I(:,:,2) = I_over_time(:,end-4000:end); 
                O(:,:,2) = ones(10000, 4001, 1) - S(:,:,2) - I(:,:,2); 
load('Gt02');   S(:,:,3) = S_over_time(:,end-4000:end); 
                I(:,:,3) = I_over_time(:,end-4000:end); 
                O(:,:,3) = ones(10000, 4001, 1) - S(:,:,3) - I(:,:,3); 
load('Gt03');   S(:,:,4) = S_over_time(:,end-4000:end); 
                I(:,:,4) = I_over_time(:,end-4000:end); 
                O(:,:,4) = ones(10000, 4001, 1) - S(:,:,4) - I(:,:,4); 
load('Gt04');   S(:,:,5) = S_over_time(:,end-4000:end); 
                I(:,:,5) = I_over_time(:,end-4000:end); 
                O(:,:,5) = ones(10000, 4001, 1) - S(:,:,5) - I(:,:,5); 
load('Gt05');   S(:,:,6) = S_over_time(:,end-4000:end); 
                I(:,:,6) = I_over_time(:,end-4000:end); 
                O(:,:,6) = ones(10000, 4001, 1) - S(:,:,6) - I(:,:,6); 
load('Gt06');   S(:,:,7) = S_over_time(:,end-4000:end); 
                I(:,:,7) = I_over_time(:,end-4000:end); 
                O(:,:,7) = ones(10000, 4001, 1) - S(:,:,7) - I(:,:,7);  
load('Gt07');   S(:,:,8) = S_over_time(:,end-4000:end); 
                I(:,:,8) = I_over_time(:,end-4000:end); 
                O(:,:,8) = ones(10000, 4001, 1) - S(:,:,8) - I(:,:,8);  
load('Gt08');   S(:,:,9) = S_over_time(:,end-4000:end); 
                I(:,:,9) = I_over_time(:,end-4000:end); 
                O(:,:,9) = ones(10000, 4001, 1) - S(:,:,9) - I(:,:,9); 
load('Gt09');   S(:,:,10) = S_over_time(:,end-4000:end); 
                I(:,:,10) = I_over_time(:,end-4000:end); 
                O(:,:,10) = ones(10000, 4001, 1) - S(:,:,10) - I(:,:,10); 
load('Gt1');    S(:,:,11) = S_over_time(:,end-4000:end); 
                I(:,:,11) = I_over_time(:,end-4000:end); 
                O(:,:,11) = ones(10000, 4001, 1) - S(:,:,11) - I(:,:,11); 

%% Find Assortativity values
for j = 1:11
    for i = 1:4001
        Matrix = [  [sum(S(:,i,j)'*(S(:,i,j)'*A)'), sum(S(:,i,j)'*(I(:,i,j)'*A)')]; ...
                    [sum(I(:,i,j)'*(S(:,i,j)'*A)'), sum(I(:,i,j)'*(I(:,i,j)'*A)')]];
        M = Matrix/(sum(Matrix, 'all'));
        R_over_global_over_time(j, i) = (trace(M(:,:)) - sum(sum(M(:,:),1)'.*sum(M(:,:),2)))/(1 - sum(sum(M(:,:),1)'.*sum(M(:,:),2)));
        II_over_global_over_time = M(2,2);
        SS_over_global_over_time = M(1,1);
    end
end

%% Figure 3b (homogeneous)

R_mean = mean(R_over_global_over_time,2);

fig1 = figure(); %relative assortativity
plot(0:0.1:1, R_mean/R_mean(end), 'LineWidth', 1.5);
hold on

% fig2 = figure(); %absolute assortativity
% plot(0:0.1:1, R_mean, 'LineWidth', 1.5);
% hold on