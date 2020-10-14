tspan = 0:100;
c_0 = 100;
weights = 1;
target = 10;
trials = 1000;
source_bound = 500;
v = 1;
m = -v.^2/2;

source_loc = [8,4]';
sinks = [7,1]';
[H,time,C,k_i,Qf,Qi,k_final,edgeidx,~,~] = hu_cai_model_with_fixed_sinks_and_sources(0,1,source_loc,sinks,source_bound,tspan,c_0,weights,m,v);

sum(k_i.^1/2)
sum(k_final.^1/2)
%generating_data_hu_cai(0,1,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m,v);
right_edge_flows = zeros(current_config,total_trials);
left_edge_flows = zeros(current_config,total_trials);
total_edgelist = (abs(total_Q_matrix) > 0)';
% figure(1)
% HWidths = max(10*k_i/max(k_i),10^-6);
% num_nodes = numnodes(H);
% p = plot(H,'Layout','force','LineWidth',HWidths);
% p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
% figure(2)
% HWidths = max(10*k_final/max(k_final),10^-6);
% num_nodes = numnodes(H);
% p = plot(H,'Layout','force','LineWidth',HWidths);
% p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
for i=1:current_config
    for j=1:total_trials
        if total_edgelist(j,:) == total_configs(i,:)
            right_edge_flows(i,j) = total_init_Q_matrix(2,j)/total_init_Q_matrix(5,j);
            left_edge_flows(i,j) = total_init_Q_matrix(10,j)/total_init_Q_matrix(5,j);
        end
    end
end
% right_edge_flows = right_edge_flows(right_edge_flows ~= 0);
% left_edge_flows = left_edge_flows(left_edge_flows ~= 0);
figure(2)
hold on
scatter(right_edge_flows(1,:), left_edge_flows(1,:),'DisplayName', 'middle_left')
scatter(right_edge_flows(2,:), left_edge_flows(2,:),'DisplayName', 'middle')
scatter(right_edge_flows(3,:), left_edge_flows(3,:),'DisplayName', 'left_right')
scatter(right_edge_flows(4,:), left_edge_flows(4,:),'DisplayName', 'middle_right')
legend;
            
        