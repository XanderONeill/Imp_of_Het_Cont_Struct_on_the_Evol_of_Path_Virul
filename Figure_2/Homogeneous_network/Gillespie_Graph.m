function [G, dt] = Gillespie_Graph(G, params)

b = params(1); d = params(2); G_T = params(3);
beta = G.Nodes.Transmission';
K = G.numnodes;

S = G.Nodes.Status == "S";      S = S';
I = G.Nodes.Status == "I";      I = I';
O = G.Nodes.Status == "O";      O = O';

%% Event probabilities
birth = b*(S*(O*G.adjacency)');
sus_death = d*sum(S);
inf_death = d*sum(I) + sum(G.Nodes.Virulence);
trans = (1 - G_T)*(beta)*(S*G.adjacency')' + G_T*sum(beta*outdegree(G))*sum(S)/K;

R_T = birth + sus_death + inf_death + trans;

%% Gillespie algorithm
dt = exprnd(1/R_T);
U = rand;
if (U <= birth/R_T)                             %birth event
    % S*A1 - for each node, how many neighbours are susceptible
    % O.*(S*A1)' - for each empty node, number of susceptible connections

    % splitting the [0 birth/R_T] interval into further categories to
    % determine which empty node now becomes susceptible
    prob_vec = [0, cumsum(b*(O.*(S*G.adjacency)))]/R_T;
%     f = find(prob_vec < U);
%     ind = f(end); % ind denotes the index of the empty node that will become sus
    f = prob_vec<U;
    ind = sum(f);
    
    G.Nodes.Status(ind) = "S";

elseif (U <= (birth+sus_death)/R_T)             %susceptible death
    % splitting the [birth/R_T (birth+sus_death)/R_T] interval into 
    % further categories to determine which susceptible node dies
    % if a node is not susceptible then the entry in the vector S will be
    % 0, and therefore cannot be chosen. if a node is susceptible then its
    % entry is 1, and is then weighted by its natural death.
    
    prob_vec = [0, cumsum(d*S)]/R_T + (birth)/R_T;
%     f = find(prob_vec < U);
%     ind = f(end);
    f = prob_vec<U;
    ind = sum(f);

    G.Nodes.Status(ind) = "O";

elseif (U <= (birth+sus_death+inf_death)/R_T)   %infected death
    % splitting the [(birth+sus_death)/R_T (birth+sus_death+inf_death)/R_T] interval into 
    % further categories to determine which infected node dies
    % if a node is not infected then its entry is 0, and therefore cannot
    % be chosen. if a node is infected then its entry is 1, and is then
    % weighted by the natural death and virulence of that individual.

    prob_vec = [0, cumsum(d*I + G.Nodes.Virulence')]/R_T + (birth + sus_death)/R_T;
%     f = find(prob_vec < U);
%     ind = f(end);
    f = prob_vec<U;
    ind = sum(f);
    
    G.Nodes.Status(ind) = "O";

    % update virulences and transmissions
    G.Nodes.Transmission(ind) = 0;
    G.Nodes.Virulence(ind) = 0;

else                                            %transmission event
    %%% ----------------------------------------
    % process for finding out which infected individual is transmitting
    % the infection

    % follow the supplementary material for details on how this works. it
    % becomes considerably easier to understand this code if you have the
    % transmission terms (from Supp Mat) in front of you. 
    
    %FD
    % LOCAL TRANS - (1-p)*(beta.*I).*(S*G.adjacency')
    % beta.*I - weight each infected individual by their respective
    % transmission coefficient
    % S*A2' - a vector where each entry is the number of susceptible
    % neighbours (irregular neighbours here) that a given node has. We take
    % the transpose of A2, rather than just A2, as we are looking for the
    % outgoing neighbours of a given node, which may be different from the
    % incoming neighbours of a given node. 

    % GLOBAL TRANS - p*sum(S)/K*sum((beta.*I)'.*G.adjacency, 2)'
    % sum(S)/K - frequency-dependence; proportion of susceptibles in the
    % grid to total grid capacity (i.e. chance an infected comes into
    % global contact with a susceptible)
    % beta.*I - again, weight each infected by their transmission coeff
    % sum((beta.*I).*A2,2)' - further weight the chance of transmission
    % based on an infected individuals degree/connectivity. The sum(:,2)
    % denotes the column summation, and the transpose is used to make sure
    % the vector has the same orientation as the vector in LOCAL TRANS

    prob_vec = [0, cumsum((1-G_T)*beta.*(S*G.adjacency') + G_T*beta.*outdegree(G)'*sum(S)/K)]/R_T + (birth + sus_death + inf_death)/R_T;

%     f = find(prob_vec < U);
%     ind = f(end); 
    f = prob_vec<U;
    ind = sum(f);

    % ind denotes the index of the infected individual that will be
    % transmitting the infection to a susceptible individual.

    %%% ----------------------------------------
    % process for finding out which susceptible individual is being
    % infected

    % again, using the supplementary material helps here. We've chosen an
    % individual, i = ind, and so now we just need to determine which
    % susceptible they infect. 

    % LOCAL TRANS - beta(ind)*((1-p)*A2(ind,:)*S'
    % I(ind) = 1, so we don't need this in the formulation
    % beta(ind) - technically don't need this in either, as beta's will
    % cancel, but leave it in for completeness
    % A2(ind,:)*S' - the susceptible neighbours of infected node I(ind).
    
    % GLOBAL TRANS - p/K*sum(S)*sum(A2(ind,:)))
    % sum(S)/K - frequency-dependence
    % sum(A2(ind,:)) - the outgoing degree of infected node I(ind) (i.e.
    % connectivity)  

    local_nodes = zeros(1, length(S));
    local_nodes(successors(G,ind)) = 1;
    local_sus = local_nodes.*S;

    sus_prob = beta(ind)*((1 - G_T)*sum(S(successors(G,ind))) + G_T/K*sum(S)*outdegree(G, ind));  
    prob_vec_2 = [0, cumsum(beta(ind)*(1 - G_T)*local_sus + beta(ind)*outdegree(G,ind)*G_T/K*S)]*(prob_vec(ind+1) - prob_vec(ind))/sus_prob + prob_vec(ind);

%     f2 = find(prob_vec_2 < U);
%     ind2 = f2(end);
    f2 = prob_vec_2<U;
    ind2 = sum(f2);

    % ind2 denotes the index of the susceptible that will be infected from
    % the m*n nodes.

    %%% -----------------------------------------
    % set the index of that susceptible to be 0 (they are no longer sus)
    % and change it to infected
    
    G.Nodes.Status(ind2) = "I";

    %%% -----------------------------------------
    % set the alpha type of that new infected to match the individual
    % that infected

    G.Nodes.Virulence(ind2) = G.Nodes.Virulence(ind);

    %%% -----------------------------------------
    % chance of evolution to occur
    
    U_E = rand;
    if (U_E<=0.005) && (G.Nodes.Virulence(ind2)>0) % 0 here is alpha_min in trade-off
        G.Nodes.Virulence(ind2) = G.Nodes.Virulence(ind2) - 0.01;
    elseif (U_E>0.995) && (G.Nodes.Virulence(ind2)<5)   % 5 here is alpha_max in trade-off
        G.Nodes.Virulence(ind2) = G.Nodes.Virulence(ind2) + 0.01;
    end

    G.Nodes.Transmission(ind2) = trade_off_function(G.Nodes.Virulence(ind2));
end  
end