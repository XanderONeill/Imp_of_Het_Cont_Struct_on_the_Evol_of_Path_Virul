function A = IrregularAdj(m,n,k,mu)
% Regular network adjacency matrix (used for birth processes)

%m - number of rows in the network
%n - number of columns in the network
%k - shape parameter of the gamma distribution
%mu - mean of the gamma distribution
theta = mu/k;

% d_vec -   produces a vector of degrees for each node. These degrees are
%           taken from a discrete gamma distribution with shape parameter, k, and
%           scale parameter, theta. Theta is always chosen such that the mean
%           connectivity for different values of shape parameter is always constant. 
%           The while loop exists as there is a possibility that a randomly sampled
%           degree is greater than the number of nodes (which would therefore mean
%           that node can't connect to all of its required neighbours). In this case,
%           we repeat the degree selection until all nodes have degree less than the
%           maximum number of neighbours
d_vec = gamrnd(k, theta, 1, m*n);
while max(d_vec)>m*n-1
    d_vec = gamrnd(k, theta, 1, m*n);
end
d_vec = gamma_rounding(d_vec, k, mu);

A1 = RegularAdj(m,n);   % A1 - adjacency matrix for the birth process
A = zeros(m*n,m*n);     % A - irregular adjacency matrix

f = find(d_vec>0);
l = length(f);
for i = 1:l
    p = floor((-2 + sqrt(4 + 8*d_vec(f(i))))/4);
    e = d_vec(f(i))- (2*p^2+2*p);

    sumMat = zeros(m*n,m*n);
    for j = 0:p
        sumMat = sumMat + A1^j;
    end
    A(f(i), sumMat(f(i),:)>0) = 1;
    A(f(i), f(i)) = 0;

    Ap1 = A1^(p+1) + sumMat;
    Ap_exclusive = setxor(find(Ap1(f(i),:)>0), find(sumMat(f(i),:)>0));
    r = randsample(Ap_exclusive, e);
    A(f(i), r) = 1;
end

A = sparse(A);
end
        





