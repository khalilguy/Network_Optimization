function [k,I,it,F,p,num_nodes,G] = create_digraph(s,t,weights,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
G = digraph(s,t,weights);
num_nodes = numnodes(G);
it = incidence(G);
I = (it)';

k = G.Edges.Weight;
F = sparse(num_nodes,1);

for i = 1: size(source_loc)
    F(source_loc(i,1),1) = source_boundary_conditions(i,1);
end


for i = 1 : size(sink_loc,1)
    F(sink_loc(i,1),1)= sink_boundary_conditions(i,1);
end

 p = randi(num_nodes);
 F(p) = 0;
end

