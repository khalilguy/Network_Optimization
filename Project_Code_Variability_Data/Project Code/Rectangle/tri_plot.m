function intersection_of_tri_configs_and_unique_config_edges_from = tri_plot(size_graph_x, size_graph_y,source_loc,sinks,source_bound,tspan,c_0,number_of_points)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sigma = 2;
mu = sigma^2/2;

load possible_configurations_hu_cai2by2.mat unique_k_configs
load particular_configuration_hu_cai_2by2.mat tri_configs
sample_k = unique_k_configs(1,:).*lognrnd(mu,sigma,1,length(unique_k_configs(1,:)));
t1 = linspace(0,1,number_of_points); t2 = linspace(0,1,number_of_points);
[x,y] = meshgrid(t1,t2);
t3 = 1-x-y;
isok = t3>=0;
t1 = x(isok); 
t2 = y(isok); 
t3=t3(isok);
new_configs = zeros(length(t3),length(sample_k));

for i = 1:length(t3)
    k_initial = t1(i,1)*tri_configs(1,:) + t2(i,1)*tri_configs(2,:) + t3(i,1)*tri_configs(3,:);
    [H,time,C1,k_i,Qf,Qi,k_final,edgeidx,Q_flow_direction,changed_flow_direction] = hu_cai_model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound,tspan,c_0,k_initial);
    new_configs(i,:) = k_final;

end

tri_configs_edges = tri_configs > source_bound/length(sample_k);
tri_configs_idx = zeros(length(new_configs),1);
config_edges = new_configs > source_bound/length(sample_k);
[unique_config_edges_from_tri_configs_optimizations] = unique(config_edges,'rows');
for i = 1:size(config_edges,1)
    for j = 1:size(unique_config_edges_from_tri_configs_optimizations,1)
        if config_edges(i,:) == unique_config_edges_from_tri_configs_optimizations(j,:)
            tri_configs_idx(i,1) = j;
        end
    end
end

h = 0;
z = zeros(size(x));
for i = 1:size(x)
    for j = 1: size(y)
        if isok(i,j) == 1
            h = h + 1;
            z(j,i) = tri_configs_idx(h,1);
        end
    end
end

intersection_of_tri_configs_and_unique_config_edges_from = intersect(tri_configs_edges,unique_config_edges_from_tri_configs_optimizations,'rows');
tri = delaunay(x,y);
trisurf(tri,x,y,z)
end

