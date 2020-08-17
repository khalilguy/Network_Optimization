function Q = calculate_flows(k,I,it,F,p,num_nodes)
%UNTITLED Summary of this function goes here
%   k is conductances on each of the edges, I is the tranpose of incidence matrix of the graph G,it
%   is original matrix given by the function incidence(G), F is a vector
%   for the sum of flows at a node(zero everywhere except for the sources
%   and the sinks), p is a random row index where we set the pressure to be
%   zero, and num_nodes is the number of nodes in the graph

K = diag(k);
R = it*K*I;


randpressure = sparse(1,num_nodes);
randpressure(1,p) = 1;
R(p,:) = randpressure; %This is done because R is singular
P = R\F;

Q = K*I*P;
end

