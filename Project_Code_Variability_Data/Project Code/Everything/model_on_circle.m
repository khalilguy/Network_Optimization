function [Q,k,H] = model_on_circle(source_loc,sinks,source_bound,radius,tol)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[out_nodes, in_nodes,sink_loc, new_nodes,num_nodes] = hexagonal_graph_fixed_sinks(10,10,source_loc,sinks);
[C, out_nodes, in_nodes, location, circle] = create_circle(new_nodes, num_nodes, radius, tol);
sink_bound = zeros(size(sink_loc,1),1);
source_bounds = zeros(size(source_loc,1),1);

for i = 1:size(source_loc,1)
    source_bounds(i,1) = source_bound/size(source_bounds,1);
end

for i = 1:size(sink_loc,1)
    sink_bound(i,1) = sum(source_bound)/size(sink_bound,1)*-1;
end

[Q,k,H] = bohn_magnasco(out_nodes,in_nodes,source_loc,sink_loc,source_bounds,sink_bound);
end

