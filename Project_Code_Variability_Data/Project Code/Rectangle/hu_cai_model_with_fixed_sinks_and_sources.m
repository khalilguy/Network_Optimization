function [H,time,C,k_i,Qf,Qi,k_final,edgeidx,Q_flow_direction,changed_flow_direction] = hu_cai_model_with_fixed_sinks_and_sources(size_graph_x, size_graph_y,source_loc,sinks,source_bound,tspan,c_0,weights,m,v)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[out_nodes, in_nodes,sink_loc] = hexagonal_graph_fixed_sinks(size_graph_x,size_graph_y,source_loc,sinks);

sink_bound = zeros(size(sink_loc,1),1);
source_bounds = zeros(size(source_loc,1),1);


for i = 1:size(source_loc,1)
    source_bounds(i,1) = source_bound/size(source_bounds,1);
end

for i = 1:size(sink_loc,1)
    sink_bound(i,1) = sum(source_bound)/size(sink_bound,1)*-1;
end


[H,time,C,k_i,Qf,Qi,k_final,Q_flow_direction,changed_flow_direction]= hu_cai(tspan,c_0,out_nodes,in_nodes ,source_loc,sink_loc,source_bounds,sink_bound,weights,m,v);
    
edgeidx = findedge(H,out_nodes, in_nodes);
end
