%%
%script used to run all the main functions. Most functions require the
%source and sink locations and change based on the dimensions of the
%network. Below are are few of the network dimensions that were used and
%there correspond source and sink locations semetrically placed at opposite
%ends of the lattice
tic
tspan = 0:100;
c_0 = 1;
weights = 1;
target = 5;
trials = 1000;
source_bound = 500;


v = 1:2:10;
m = -v.^2/2;

%0,1

for i = 1:length(v)
    source_loc = [8,4]';
    sinks = [7,1]';
    generating_data_hu_cai(0,1,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m(i),v(i));
    analysis_hu_cai(0,1,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m(i),v(i));
    close all
end


% for i = 1:length(v)
%     source_loc = [32 36 40]';
%     sinks = [1 7 11]';
%     possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m(i),v(i));
%     particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m(i),v(i));
%     close all
% end 

% 0,1
% tic
% source_loc = [8,4]';
% sinks = [7,1]';
% possible_configurations_hu_cai(0,1,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(0,1,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% first_part_time = toc;
% 
% % 2,2
% tic
% source_loc = [32 40]';
% sinks = [1 11]';
% possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% second_part_time = toc;
% 
% tic;
% source_loc = [32 36 40]';
% sinks = [1 7 11]';
% possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% third_part_time = toc;
% 
% tic;
% source_loc = [32 36 40]';
% sinks = [1 11]';
% possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% fourth_part_time = toc;
% % 
% tic;
% source_loc = [32 36 40]';
% sinks = 7;
% possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% fifth_part_time = toc;
% % 
% tic
% source_loc = [32 40]';
% sinks = [1 7 11]';
% possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% sixth_part_time = toc;
% % 
% tic
% source_loc = 36;
% sinks = [1 7 11]';
% possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% seventh_part_time = toc;
% 
% 
% %2,3
% tic
% source_loc = [32 36 40]';
% sinks = [1 7 11]';
% possible_configurations_hu_cai(2,3,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,3,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% eigth_part_time = toc;
% % % source_loc = [15,11,7,1]';
% % % sinks = [52,48,44,40]';
% % 
% % %2,5
% % % source_loc = 14;
% % % sinks = 65;
% %2,9
% tic
% source_loc = 147;
% sinks = 42;
% possible_configurations_hu_cai(2,9,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,9,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% close all
% nineth_part_time = toc;
% entire_time = toc;
% % % source_loc = [27,23,19,15]';
% % % sinks = [112,108,104,100]';
% % 
% save('times.dat', 'first_part_time', 'second_part_time', 'third_part_time', 'fourth_part_time', 'fifth_part_time', 'sixth_part_time', 'seventh_part_time', 'eighth_part_time', 'nineth_part_time'); 
% % % 
% [total_configs,config_nums,total_trials,absolute_A,total_frequency] = possible_configurations(2,2,source_loc,sinks,source_bound,trials,target);
% particular_configuration(2,2,source_loc,sinks,source_bound,trials,target);
% [A2,f3,f4] = energy_and_frequency_vs_edges(2,2,source_loc,sinks,source_bound,trials);


%%
%Looking at initial conductances
% trials = 1000;
% [final_configs, initial_distance, I] = initial_k(2,2,source_loc,sinks,source_bound,trials);
% [Q,k,H,changed_direction,k_i,edgeidx,Qf] = model_with_fixed_sinks_and_sources(2,2,source_loc,sinks,source_bound);
% widths_i = 10*k_i/max(k_i);
% widths = 10*k/max(k);
% num_nodes = numnodes(H);
% figure
% subplot(1,2,1)
% p = plot(H,'Layout','force');
% p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
% subplot(1,2,2)
% p = plot(H,'Layout','force','LineWidth',widths);
% p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
%conductances which are zero correspond to flows which are zero, thus a
%unique configuration
 
%%
% [out_nodes, in_nodes,sink_loc] = hexagonal_graph_fixed_sinks(3,9,source_loc,sinks);
% tspan = 0:40;
% c_0 = 1;
% num_edges = length(out_nodes);
% weights = ones(1,num_edges) + .5*rand(1,num_edges);
% [H,~,~,~,Q,~,k_final,~,~,~] = hu_cai_model_with_fixed_sinks_and_sources(3,9,source_loc,sinks,source_bound,tspan,c_0,weights);
% for i = 1:trials
% [k,H,time,C,k_i,Qf,Qi,k_final,edgeidx] = hu_cai_model_with_fixed_sinks_and_sources(2,2,source_loc,sinks,source_bound,tspan,c_0,weights);
% % C = real(C);
% figure(1)
% plot(time, C(:,:))
% figure(2)
% HWidths = max(5*k_final/max(k_final),10^-6);
% plot(H,'Layout','force','EdgeLabel',k_final,'LineWidth',HWidths)
% end
% 
% figure(3)
% HWidths = max(5*k_final/max(k_final),10^-6);
% num_nodes = numnodes(H);
% p = plot(H,'Layout','force','LineWidth',HWidths);
% p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);

%%
tspan = 0:100;
c_0 = 1;
weights = 1;
% % [A] = hu_cai_energy_and_frequency_vs_edges(2,2,source_loc,sinks,source_bound,trials,tspan,c_0,weights);
%[current_config,configs_num,total_configs,total_trials,absolute_A,total_frequency,total_Q_matrix,total_energy_array] = possible_configurations_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);
% particular_configuration_hu_cai(2,2,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights);

%%

% tspan = 0:100;
% c_0 = 1;
% tri_plot(2,2,source_loc,sinks,source_bound,tspan,c_0,100)
%%
% sigma = 2;
% mu = sigma^2/2;
% load possible_configurations_hu_cai2by2.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges unique_k_configs unique_Q_configs
% load('particular_configuration_hu_cai_2by2.mat')
% sample_k = unique_k_configs(1,:).*lognrnd(mu,sigma,1,length(unique_k_configs(1,:)));
% tspan = 0:100;
% c_0 = 1;
% 
% 
% t2 = linspace(0,1,100);
% t1 = 1 - t2;
% 
% pw_size = size(top_unique_k_configurations,1);
% pairwise_winners = triu(ones(pw_size,pw_size));
% pairwise_winning_scores = triu(ones(pw_size,pw_size));
% 
% intersection_of_top_configs_with_config_edges = triu(ones(pw_size,pw_size));
% top_unique_k_config_edges = top_unique_k_configurations > source_bound/length(sample_k);
% 
% for i = 1 : size(top_unique_k_configurations,1)
%     for j = 1 : size(top_unique_k_configurations,1)
%         if i ~= j
%             if pairwise_winners(i,j) ~= 0
%                 new_configs = zeros(length(t1),length(sample_k));
%                 for k = 1:length(t1)
%                     k_initial = t1(1,k)*top_unique_k_configurations(i,:) + t2(1,k)*top_unique_k_configurations(j,:);
%                     [H,time,C1,k_i,Qf,Qi,k_final,edgeidx,Q_flow_direction,changed_flow_direction] = hu_cai_model_with_fixed_sinks_and_sources(2,2,source_loc,sinks,source_bound,tspan,c_0,k_initial);
%                     new_configs(k,:) = k_final;
%                 end
%                 
%                 config_edges = new_configs > source_bound/length(sample_k);
%                 intersection_new_configs_top_unique_k_configs = intersect(config_edges,top_unique_k_config_edges,'rows');
%                 [unique_config_edges_from_tri_configs_optimizations] = unique(config_edges,'rows');
%                 pair_configs_idx = zeros(size(config_edges,1),1);
% 
%                 for l = 1:size(config_edges,1)
%                     for m = 1:size(unique_config_edges_from_tri_configs_optimizations,1)
%                         if config_edges(l,:) == unique_config_edges_from_tri_configs_optimizations(m,:)
%                             pair_configs_idx(l,1) = m;
%                         end
%                     end
%                 end
% 
%                 end_element = pair_configs_idx(1,1);
%                 flipped_idx = flip(pair_configs_idx);
%                 other_end_element = flipped_idx(1,1);
%                 end_length = 0;
%                 other_end_length = 0;
% 
%                 for q = 1:length(pair_configs_idx)
%                     if pair_configs_idx(q,1) ~= end_element
%                         break
%                     end
%                     if pair_configs_idx(q,1) == end_element
%                         end_length = end_length + 1;
%                     end
%                 end
%                 for q = 1:length(pair_configs_idx)
%                      if flipped_idx(q,1) ~= other_end_element
%                         break
%                     end
%                     if flipped_idx(q,1) == other_end_element
%                         other_end_length = other_end_length + 1;
%                     end
% 
%                 end
%                 winning_score = max(end_length, other_end_length);
%                 
%                  if end_length == other_end_length
%                     pairwise_winners(i,j) = 0;
%                     pairwise_winning_scores(i,j) = 0;
%                  elseif end_length ==  winning_score
%                     pairwise_winners(i,j) = i;
%                     pairwise_winning_scores(i,j) = end_length;
%                  else
%                     pairwise_winners(i,j) = j;
%                     pairwise_winning_scores(i,j) = other_end_length;
%                 end
%  
%             end  
%         end
%     end
% end
%%
% load possible_configurations_hu_cai2by2.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges unique_k_configs unique_Q_configs
% load('particular_configuration_hu_cai_2by2.mat')

% path4 = 'hu_cai_particular_plots\configurations_graphs';
% weights = 1;
% [H,~,~,~,~,~,~,~,~,~] = hu_cai_model_with_fixed_sinks_and_sources(2,2,source_loc,sinks,source_bound,tspan,c_0,weights);
% 
% 
% for i = 1:23
%     figure
%     HWidths = 10*ordered_unique_k_configs(i,:)/max(ordered_unique_k_configs(i,:));
%     num_nodes = numnodes(H);
%     p = plot(H,'Layout','force','LineWidth',HWidths);
%     p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
%     saveas(gca, fullfile(path4,['config' int2str(i) '.fig']), 'fig');
%     saveas(gca, fullfile(path4,['config' int2str(i) '.jpeg']), 'jpeg')
% 
% end

