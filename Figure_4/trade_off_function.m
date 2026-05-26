function f = trade_off_function(alpha)      % this seems self-explanatory
a = -0.75;
beta_min = 1;
beta_max = 6;
alpha_min = 0;
alpha_max = 5;
%f = 3*alpha;
f = beta_min + ((beta_max - beta_min)*(1 - (alpha_max - alpha)/(alpha_max - alpha_min)))./(1 + a*(alpha_max - alpha)/(alpha_max - alpha_min));
end
