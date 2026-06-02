function indices = add_new_connections(G_r, ind, k, mu)
% function to find new indices that the new susceptiple connects to

d = gamrnd(k, mu/k, 1, 1);
while d > length(G_r.adjacency)-1
    d = gamrnd(k, mu/k, 1, 1);
end
d = round(d);

p = floor((-2 + sqrt(4 + 8*d))/4);
e = d - (2*p^2+2*p);

indices = nearest(G_r, ind, p);
[nn, dist] = nearest(G_r, ind, p+1);
indices = [indices; randsample(nn(dist==p+1), e)];

end