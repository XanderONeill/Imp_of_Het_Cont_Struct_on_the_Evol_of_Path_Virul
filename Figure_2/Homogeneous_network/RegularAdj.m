function A = RegularAdj(m,n)
% Regular network adjacency matrix (used for birth processes)

%m - number of rows in the network
%n - number of columns in the network

A = zeros(m*n,m*n);    % adjacency matrix for the birth process

for i = 1:m*n               % go through each node in the network
    if (mod(i,n) ~= 0)      % if the current node isn't in the final column, then make a connection to the adjacent (right) neighbour
        A(i,i+1) = 1;
    else                    % if the current node is in the final column, then make a connection to the node in the first column (of that row)
        A(i,i-n+1) = 1;
    end

    if (i <= (m-1)*n)       % if the current node isn't in the final row, then make a connection to the adjacent (down) neighbour
        A(i,i+n) = 1;
    else
        A(i,i-(m-1)*n) = 1; % if the current node is in the final row, then make a connection to the node in the first row (of that column)
    end
end

A = A + A';         % the regular network is symmetrical (i.e. connections
                    % right and down are also connections left and up)
                    % and so add the transpose to the original adjacency
                    % matrix to get all left and up connections
A = sparse(A);
end

