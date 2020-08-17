function [Q,k,H] = model_with_uneven_sinks_and_sources(size_graph_x, size_graph_y,source_loc,sink_loc,source_bounds,sink_bounds)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[out_nodes, in_nodes,sink_loc] = hexagonal_graph_fixed_sinks(size_graph_x,size_graph_y,source_loc,sink_loc);

[Q,k,H] = bohn_magnasco(out_nodes,in_nodes,source_loc,sink_loc,source_bounds,sink_bounds);
end



