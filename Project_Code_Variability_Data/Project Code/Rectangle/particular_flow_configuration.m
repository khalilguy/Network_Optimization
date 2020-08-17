function [flow_configuration_frequency,total_flow_configs] = particular_flow_configuration(size_graph_x, size_graph_y,source_loc,sinks,source_bound,trials,target)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% load possible_configurations0by1.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges changed_direction_array total_flow_configs
load possible_configurations2by2.mat total_configs total_flow_configs total_flow_trials total_Qf_matrix changed_direction_array absolute_A
% load possible_configurations2by3.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges changed_direction_array total_flow_configs

flow_configuration_frequency = zeros(size(total_flow_configs,1),1);
ordered_flow_configs = zeros(size(total_flow_configs,1), size(total_flow_configs,2));
current_config = size(total_flow_configs,1);

m = 0;
for j = 1:size(absolute_A ,1)
    for i = 1:size(total_flow_configs,1)
        if absolute_A(j,1) == sum(total_configs(i,:))
            m = m + 1;
            ordered_flow_configs(m,:) = total_flow_configs(i,:);
        end
    end
end


flow_distance_array = zeros(current_config*(current_config -1),1);
total_flow_configs_idx = 1:current_config;

flow_distance_idx_array = zeros(current_config*(current_config -1), 1);
flow_distance_matrix = zeros(current_config,current_config);

avg_flow_difference = zeros(current_config,1);
k = 1;



for i = 1:current_config 
    config_i = total_flow_configs(i,:);
    for  p = 1:current_config
        if p ~= i
            config_j = total_flow_configs(p,:);
            different_flow = 0;
            for j = 1:length(config_i)   
                if config_i(1,j) ~= config_j(1,j)
                    avg_flow_difference(i,1) = avg_flow_difference(i,1) + 1;
                    different_flow = different_flow + 1;
                end
            end
            flow_distance_matrix(i,p) = different_flow;
            flow_distance_idx_array(k,1) = i;
            flow_distance_array(k,1) = different_flow;
            k = k + 1;
        end
    end
    avg_flow_difference(i,1) = avg_flow_difference(i,1)/(current_config - 1);
end


%%
% min_flow_distance = zeros(current_config,1);
% max_flow_distance = zeros(current_config,1);
% min_flow_idx  = 0;
% max_flow_idx = 0;
% 
% 
% l =1;
% for i = 1:current_config -1:size(distance_array,1)
%     
%     flow_dist_config = flow_distance_array(i: i + current_config - 2,1);
%     
%     [min_flow_distance(l,1), ~] = min(flow_dist_config,[],1);
%     [max_flow_distance(l,1),~] = max(flow_dist_config,[],1);
%     
%     for j = 1: length(dist_config)
%         if flow_dist_config(j,1) == min_flow_distance(l,1)
%             if min_flow_idx == 0
%                 if j>=l
%                     min_flow_idx = j+1;
%                     config_idx_flow_min = l;
%                     min_flow_configs = total_flow_configs(j+1,:);
%                 else
%                     min_flow_idx = j;
%                     config_idx_flow_min = l;
%                     min_flow_configs = total_flow_configs(j,:);
%                 end
%                     
%             elseif min_flow_idx ~= 0
%                 if j>=l
%                     min_flow_idx = [min_flow_idx j+1];
%                     config_idx_flow_min = [config_idx_flow_min l];
%                     min_flow_configs = [min_flow_configs total_flow_configs(j+1,:)];
%                 else
%                     min_flow_idx = [min_flow_idx j];
%                     config_idx_flow_min = [config_idx_flow_min l];
%                     min_flow_configs = [min_flow_configs total_flow_configs(j,:)];
%                 end
%                     
%             end
%         end
%         
%          if flow_dist_config(j,1) == max_flow_distance(l,1)
%             if max_flow_idx == 0
%                 if j>=l
%                     max_flow_idx = j+1;
%                     config_idx_flow_max = l;
%                     max_flow_configs = total_flow_configs(j+1,:);
%                 else
%                     max_flow_idx = j;
%                     config_idx_flow_min = l;
%                     max_flow_configs = total_flow_configs(j,:);
%                 end
%                     
%             elseif max_flow_idx ~= 0
%                 if j>=l
%                     max_flow_idx = [max_flow_idx j+1];
%                     config_idx_flow_max = [config_idx_flow_max l];
%                     max_flow_configs = [max_flow_configs total_flow_configs(j,:)];
%                 else
%                     max_flow_idx = [max_flow_idx j];
%                     config_idx_flow_max = [config_idx_flow_max l];
%                     max_flow_configs = [max_flow_configs total_flow_configs(j,:)];
%                 end
%                     
%             end
%          end
%     end
%     
%     l = l+1;
% end
% 
% 
% min_flow_freq = zeros(current_config,1);
% max_flow_freq = zeros(current_config,1);
% 
% 
% for i = 1:current_config
%     
%     for j = 1:length(min_flow_idx)
%         if total_flow_configs_idx(1,i) == min_flow_idx(1,j)
%             min_flow_freq(i,1) = min_flow_freq(i,1) + 1;
%         end
%     end
%         
%     for j = 1:length(max_flow_idx)
%         if total_flow_configs_idx(1,i) == max_flow_idx(1,j)
%             max_flow_freq(i,1) = max_flow_freq(i,1) + 1;
%         end
%     end
% end
% 
% 
% changed_direction_freq = [(total_trials - sum(changed_direction_array)) sum(changed_direction_array)];


save particular_flow_configuration2by2.mat total_flow_configs ordered_flow_configs
%%

% top_freqs = 4;
% [max_freqs, max_freqs_idx] = maxk(flow_configuration_frequency,top_freqs);
% [max_dists, max_dists_idx] = maxk(avg_distance_array,top_freqs);
% top_configurations = zeros(top_freqs,size(total_configs,2));
% 
% 
% for i = 1:top_freqs
%     top_configurations(i,:) = total_flow_configs(max_freqs_idx(i,1),:);
% end
% 
% h = 0;
% 
% seen_configs = zeros(top_freqs,size(total_configs,2));
% path3 = 'particular_plots\top_configs2by2';
% for i = 1:top_freqs
%     h = 0;
%     while h == 0
%     [~,k,H,~,~,~,Qf] = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
%         if Qf == top_configurations(i,:)
%             if ismember(c,seen_configs,'rows') == 0
%                 h = h + 1;
%                 seen_configs(i,:) = Qf';
%                 figure
%                 HWidths = 10*k/max(k);
%                 num_nodes = numnodes(H);
%                 p = plot(H,'Layout','force','LineWidth',HWidths);
%                 p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
%                 saveas(gca, fullfile(path3,['config' int2str(i) '.fig']), 'fig');
%                 saveas(gca, fullfile(path3,['config' int2str(i) '.jpeg']), 'jpeg')
%                 
%             end
%         end
%     end
% end
% 
% path1 = 'particular_plots\by_geometry';
% path2 = 'particular_plots\by_flow';
% 
% set(0,'DefaultFigureVisible','on');
% 
% 
% figure('Name','Configuration Frequency');
% plot(total_flow_configs_idx, flow_configuration_frequency)
% title('Total Frequency vs. Configuration')
% xlabel('Configuration')
% ylabel('Total Frequency')
%     
% saveas(gca, fullfile(path1,'total_cofiguration_freq.fig'), 'fig');
% saveas(gca, fullfile(path1,'total_cofiguration_freq.jpeg'), 'jpeg')
% 
% 
% figure('Name', 'Flow Difference vs Configuration');
% scatter(total_flow_configs_idx, avg_flow_difference)
% title('Flow Difference vs Configuration')
% xlabel('Configuration')
% ylabel('Flow Difference')
% saveas(gca, fullfile(path2,'avg_flow_dist'), 'fig');
% saveas(gca, fullfile(path2,'avg_flow_dist'), 'jpeg');
% 
% figure('Name', 'Flow Difference vs Configuration');
% scatter(flow_distance_idx_array, flow_distance_array)
% title('Flow Difference vs Configuration')
% xlabel('Configuration')
% ylabel('Flow Difference')
% saveas(gca, fullfile(path2,'avg_flow_dist'), 'fig');
% saveas(gca, fullfile(path2,'avg_flow_dist'), 'jpeg');
% 
% figure('Name', 'Changed Direction');
% plot(0:1, changed_direction_freq)
% title('Frequency vs Changed_Direction')
% xlabel('Changed Direction')
% ylabel('Frequency')
% saveas(gca, fullfile(path2,'changed_direction'), 'fig');
% saveas(gca, fullfile(path2,'changed_direction'), 'jpeg');
% 
% fg7 = figure('Name', 'Minimum Distance from Optimum vs Configuration');
% scatter(total_flow_configs_idx, min_flow_distance)
% title('Minimum Distance from Optimum vs Configuration')
% xlabel('Configuration')
% ylabel('Distance')
% saveas(gca, fullfile(path2,'min_flow_dist'), 'fig');
% saveas(gca, fullfile(path2,'min_flow_dist'), 'jpeg');
% 
% fg8 = figure('Name', 'Closest Configuration vs Configuration');
% scatter(config_idx_flow_min, min_flow_idx)
% title('Closest Configuration vs Configuration')
% xlabel('Configuration')
% ylabel('Closest Configuration')
% saveas(gca, fullfile(path2,'closest_config_flow'), 'fig');
% saveas(gca, fullfile(path2,'closest_config_flow'), 'jpeg');
% 
% 
% fg9 = figure('Name', 'Frequency vs. Closest Configurations');
% scatter(total_flow_configs_idx, min_flow_freq)
% title('Frequency vs. Closest Configurations')
% xlabel('Closest Configuration')
% ylabel('Frquency')
% saveas(gca, fullfile(path2,'closest_flow_freq'), 'fig');
% saveas(gca, fullfile(path2,'closest_flow_freq'), 'jpeg');
% 
% fg10 = figure('Name', 'Histogram of Frequency of Minimum Distance Frequencies');
%  [~,edges] = histcounts(log10(min_flow_freq));
% histogram(min_flow_freq,10.^edges)
% set(gca, 'xscale','log')
% title('Number of Different Configurations vs. Frequency of Minimum Distance Frequencies')
% xlabel('Frequency of Minimum Distance Frequencies')
% ylabel('Number of Different Configurations')
% saveas(gca, fullfile(path2,'closest_flow_histogram'), 'fig');
% saveas(gca, fullfile(path2,'closest_flow_histogram'), 'jpeg');
% 
% 
% fg7 = figure('Name', 'Maximum Distance from Optimum vs Configuration');
% scatter(total_flow_configs_idx, max_flow_distance)
% title('Maximum Distance from Optimum vs Configuration')
% xlabel('Configuration')
% ylabel('Distance')
% saveas(gca, fullfile(path2,'max_flow_dist'), 'fig');
% saveas(gca, fullfile(path2,'max_flow_dist'), 'jpeg');
% 
% fg8 = figure('Name', 'Farthest Configuration vs Configuration');
% scatter(config_idx_flow_max, max_flow_idx)
% title('Farthest Configuration vs Configuration')
% xlabel('Configuration')
% ylabel('Farthest Configuration')
% saveas(gca, fullfile(path2,'farthest_flow'), 'fig');
% saveas(gca, fullfile(path2,'farthest_flow'), 'jpeg');
% 
% figure('Name', 'Frequency vs. Farthest Configurations');
% scatter(total_flow_configs_idx, max_flow_freq)
% title('Frequency vs. Farthest Configurations')
% xlabel('Farthest Configuration')
% ylabel('Frquency')
% saveas(gca, fullfile(path2,'farthest_flow_freq'), 'fig');
% saveas(gca, fullfile(path2,'farthest_flow_freq'), 'jpeg');
% 
% figure('Name', 'Histogram of Frequency of Maximum Distance Frequencies');
%  [~,edges] = histcounts(log10(max_flow_freq));
% histogram(max_flow_freq,10.^edges)
% set(gca, 'xscale','log')
% title('Number of Different Configurations vs. Frequency of Maximum Distance Frequencies')
% xlabel('Frequency of Maximum Distance Frequencies')
% ylabel('Number of Different Configurations')
% saveas(gca, fullfile(path2,'farthest_flow_histogram'), 'fig');
% saveas(gca, fullfile(path2,'farthest_flow_histogram'), 'jpeg');
% 
% fg6 = figure('Name', 'Number of Edges Used vs Configuration');
% plot(total_flow_configs_idx,flow_edges_array)

end