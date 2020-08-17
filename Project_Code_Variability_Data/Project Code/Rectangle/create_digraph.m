function [k,I,it,F,p,num_nodes,G] = create_digraph(s,t,weights,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions)
%create_digraph - creates a direct graph representation of a hexagonal
%network
%s-is an array with node numbers that have flow coming out of them

%t-is an array with node numbers that have flow coming in to them

%source_loc-is an array or integer that tells the node number the flow will
%be coming out of

%sink_loc-is an array or integer that tells the node number the flow will
%be going into

%source_boundary_conditions-is an array will values that tell the amount of
%flow to be allocated to all the source locations from source_loc

%sink_boundary_conditions - is an array with values that tell the amount of
%flow the will be leaving through all the sink locations from sink_loc

%This function returns:
% G - is a directed graph created from s and t
%it - the  transpose of incidence matrix for that we want to use in the calculation of
%the flows
%I - the properly oriented incidence matrix
%num_nodes - the numnber of nodes in the graph G
%F - vector that will hold the flows for each edges in the graph G
%p - a random integer which tells the node that will be given a random
%pressure
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

