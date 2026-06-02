function A = IrregularAdj(G_R, k,mu)
% Regular network adjacency matrix (used for birth processes)

%m - number of rows in the network
%n - number of columns in the network
%k - shape parameter of the gamma distribution
%mu - mean of the gamma distribution
theta = mu/k;
m = numnodes(G_R);

% d_vec -   produces a vector of degrees for each node. These degrees are
%           taken from a discrete gamma distribution with shape parameter, k, and
%           scale parameter, theta. Theta is always chosen such that the mean
%           connectivity for different values of shape parameter is always constant. 
%           The while loop exists as there is a possibility that a randomly sampled
%           degree is greater than the number of nodes (which would therefore mean
%           that node can't connect to all of its required neighbours). In this case,
%           we repeat the degree selection until all nodes have degree less than the
%           maximum number of neighbours
d_vec = gamrnd(k, theta, 1, m);
while max(d_vec)>m-1
    d_vec = gamrnd(k, theta, 1, m);
end
d_vec = round(d_vec);

A = zeros(m);     % A - irregular adjacency matrix

for i = 1:length(d_vec)
    p = floor((-2 + sqrt(4 + 8*d_vec(i)))/4);
    e = d_vec(i) - (2*p^2+2*p);

    indices = nearest(G_R, i, p);
    A(i, indices) = 1;

    [nn, dist] = nearest(G_R, i, p+1);
    r = randsample(nn(dist==p+1), e);
    A(i, r) = 1;    
end

A = sparse(A);
end
        





