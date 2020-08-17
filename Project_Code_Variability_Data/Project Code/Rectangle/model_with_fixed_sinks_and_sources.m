function [Q,k,H,changed_direction,k_i,edgeidx,Qf,Qi,k_final,C,time] = model_with_fixed_sinks_and_sources(size_graph_x, size_graph_y,source_loc,sinks,source_bound,tspan,c_0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[out_nodes, in_nodes,sink_loc] = hexagonal_graph_fixed_sinks(size_graph_x,size_graph_y,source_loc,sinks);

sink_bound = zeros(size(sink_loc,1),1);
source_bounds = zeros(size(source_loc,1),1);
changed_direction = 0;

for i = 1:size(source_loc,1)
    source_bounds(i,1) = source_bound/size(source_bounds,1);
end

for i = 1:size(sink_loc,1)
    sink_bound(i,1) = sum(source_bound)/size(sink_bound,1)*-1;
end

type = 0;

if type == 1
    [Q,k,H,changed_direction,k_i,Qf,Qi] = bohn_magnasco(out_nodes,in_nodes,source_loc,sink_loc,source_bounds,sink_bound);
elseif type == 2
    [k,H,time,C,k_i,Qf,Qi,k_final]= hu_cai(tspan,c_0,out_nodes,in_nodes ,source_loc,sink_loc,source_bounds,sink_bound);
end
    
edgeidx = findedge(H,out_nodes, in_nodes);
end

