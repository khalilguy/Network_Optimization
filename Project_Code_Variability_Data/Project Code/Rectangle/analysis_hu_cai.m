function analysis_hu_cai(size_graph_x, size_graph_y,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m,v)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%File path when working on Windows
%filename = ['data\gd_hu_cai' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v.mat'];

%File path when working on MacOS
filename = ['data/gd_hu_cai' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v.mat'];


% load possible_configurations_hu_cai2by2.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges unique_k_configs unique_Q_configs unique_energy
load(filename, 'current_config', 'configs_num', 'total_configs', 'total_trials', 'absolute_A', 'total_frequency', 'total_Q_matrix', 'total_energy_array', 'total_used_edges', 'unique_k_configs', 'unique_Q_configs', 'unique_energy')
% load possible_configurations2by3.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges changed_direction_array total_flow_configs


total_edgelist = (abs(total_Q_matrix) > 0)';
ordered_configs = zeros(size(total_configs,1), size(total_configs,2));
ordered_unique_k_configs = zeros(size(total_configs,1), size(total_configs,2));
ordered_unique_Q_configs = zeros(size(total_configs,1), size(total_configs,2));
configuration_frequency = zeros(size(total_configs,1),1);
energy_array_idx = zeros(total_trials,1);
avg_energy_array = zeros(size(total_configs,1), 1);
ordered_energy = zeros(size(total_configs,1),1);


b = 0;
for j = 1:size(absolute_A ,1)
    for i = 1:size(total_configs,1)
        if absolute_A(j,1) == sum(total_configs(i,:))
            b = b + 1;
            ordered_unique_k_configs(b,:) = unique_k_configs(i,:);
            ordered_unique_Q_configs(b,:) = unique_Q_configs(i,:);
            ordered_configs(b,:) = total_configs(i,:);
            ordered_energy(b,1) = unique_energy(i,1);
        end
    end
end


for i = 1:size(ordered_configs,1)
    for j = 1:total_trials
        if ordered_configs(i,:) == total_edgelist(j,:)
            avg_energy_array(i,1) = avg_energy_array(i,1) + total_energy_array(j,1);
            energy_array_idx(j,1) = i;
            configuration_frequency(i,1) = configuration_frequency(i,1) + 1;
        end
    end
    avg_energy_array(i,1) = avg_energy_array(i,1)/configuration_frequency(i,1);
end


% gamma = .5;
% for i = 1:length(unique_k_configs)
%     E = 0;
%     for j = 1:size(unique_k_configs, 2)
%         E = E + (ordered_unique_Q_configs(i,j)^2/ordered_unique_k_configs(i,j) + (ordered_unique_k_configs(i,j)^gamma)*c_0);
%     end
%     ordered_energy(i,1) = E;
% end

big_E = round(ordered_energy, 4);
little_E = unique(big_E);

bins = 10;
[binned_E E_bins] = discretize(big_E,bins);
avg_freq_for_energy_level = zeros(length(little_E),1);

avg_of_binned_energy = zeros(bins,1);
avg_of_binned_configuration_freq = zeros(bins,1);

for i = 1:bins
    b = 0;
    for j = 1:length(binned_E)
        if binned_E(j,1) == i
            b = b + 1;
            avg_of_binned_energy(i,1) = avg_of_binned_energy(i,1) + big_E(j,1);
            avg_of_binned_configuration_freq(i,1) = avg_of_binned_configuration_freq(i,1) + configuration_frequency(j,1);
        end      
    end
    avg_of_binned_energy(i,1) = avg_of_binned_energy(i,1)/b;
    avg_of_binned_configuration_freq(i,1) = avg_of_binned_configuration_freq(i,1)/b;
end

% avg_freq_for_energy_level = zeros(length(little_E),1);
% 
% for i = 1:length(little_E)
%     m = 0;
%     for j = 1:length(big_E)
%         if little_E(i,1) == big_E(j,1)
%             avg_freq_for_energy_level(i,1) = avg_freq_for_energy_level(i,1) + configuration_frequency(j,1);
%             m = m + 1;
%         end
%     end
%     avg_freq_for_energy_level(i,1) = avg_freq_for_energy_level(i,1);
% end

% [ordered_freqs, ordered_freq_idx] = maxk(configuration_frequency, length(configuration_frequency));
% 
% 
% ordered_freq_configs = zeros(size(ordered_configs,1),size(ordered_configs,2));
% 
% for i = 1:length(configuration_frequency)
%     ordered_freq_configs(i,:) = ordered_configs(ordered_freq_idx(i,1),:);
% end
% 
% figure('Name','Frequency from largest to smallest')
% plot(1:size(ordered_configs,1),ordered_freqs) 


ordered_configs_idx = 1:size(ordered_configs,1);
% edges_array = zeros(current_config,1);
% 
% distance_idx_array = zeros(current_config*(current_config -1), 1);
% distance_array = zeros(current_config*(current_config-1), 1);
% distance_matrix = zeros(current_config,current_config);
% 
% avg_distance_array  = zeros(current_config,1);
% avg_energy_difference = zeros(current_config,1);
% 
% k = 1;
% 
% 
% for i = 1:current_config 
%     config_i = ordered_configs(i,:);
%     for  p = 1:current_config
%         avg_energy_difference(i,1) = avg_energy_difference(i,1) + sqrt((avg_energy_array(i,1) - avg_energy_array(p,1))^2);
%         if p ~= i
%             config_j = ordered_configs(p,:);
%             different = 0;
%             for j = 1:size(config_i,2)
%                 if config_i(1,j) ~= config_j(1,j)
%                     avg_distance_array(i,1) = avg_distance_array(i,1) + 1;
%                     different = different + 1;
%                 end
%                 
%             end
%             distance_matrix(i,p) = different;
%             distance_idx_array(k,1) = i;
%             distance_array(k,1) = different;
%             k = k + 1;
%         end
%     end
%     avg_distance_array(i,1) = avg_distance_array(i,1)/(current_config -1);
%     avg_energy_difference(i,1) = avg_energy_difference(i,1)/(current_config - 1);
%     edges_array(i,1) = sum(config_i);
% end
% 
% 
% %%
% min_distance = zeros(current_config,1);
% max_distance  = zeros(current_config,1);
% min_idx  = 0;
% max_idx = 0;
% 
% 
% 
% l =1;
% for i = 1:current_config -1:size(distance_array,1)
% 
%     dist_config = distance_array(i: i + current_config - 2,1);
%     
%     
%     [min_distance(l,1), ~] = min(dist_config,[],1);
%     [max_distance(l,1),~] = max(dist_config,[],1);
%     
%     
%     for j = 1: length(dist_config)
%         if dist_config(j,1) == min_distance(l,1)
%             if min_idx == 0
%                 if j>=l
%                     min_idx = j+1;
%                     config_idx_min = l;
%                     min_configs = ordered_configs(j+1,:);
%                 else
%                     min_idx = j;
%                     config_idx_min = l;
%                     min_configs = ordered_configs(j,:);
%                 end
%                     
%             elseif min_idx ~= 0
%                 if j>=l
%                     min_idx = [min_idx j+1];
%                     config_idx_min = [config_idx_min l];
%                     min_configs = [min_configs ordered_configs(j+1,:)];
%                     
%                 else
%                     min_idx = [min_idx j];
%                     config_idx_min = [config_idx_min l];
%                     min_configs = [min_configs ordered_configs(j,:)];
%                 end
%                     
%             end
%         end
%         
%         if dist_config(j,1) == max_distance(l,1)
%             if max_idx == 0
%                 if j>=l
%                     max_idx = j+1;
%                     config_idx_max = l;
%                     max_configs = ordered_configs(j+1,:);
%                 else
%                     max_idx = j;
%                     config_idx_max = l;
%                     max_configs = ordered_configs(j,:);
%                 end
%                     
%             elseif min_idx ~= 0
%                 if j>=l
%                     max_idx = [max_idx j+1];
%                     config_idx_max = [config_idx_max l];
%                     max_configs = [max_configs ordered_configs(j+1,:)];
%                 else
%                     max_idx = [max_idx j];
%                     config_idx_max = [config_idx_max l];
%                     max_configs = [max_configs ordered_configs(j,:)];
%                 end
%                     
%             end
%         end
%         
%     end
%     
%     l = l+1;
% end
% 
% 
% min_freq = zeros(current_config,1);
% max_freq = zeros(current_config,1);
% 
% 
% 
% for i = 1:current_config
%     for j = 1:length(min_idx)
%         if ordered_configs_idx(1,i) == min_idx(1,j)
%             min_freq(i,1) = min_freq(i,1) + 1;
%         end
%     end
%     
%     for j = 1:length(max_idx)
%         if ordered_configs_idx(1,i) == max_idx(1,j)
%             max_freq(i,1) = max_freq(i,1) + 1;
%         end
%     end
% end


%%
set(0,'DefaultFigureVisible','on');
%File path when working on Windows
% path = ['hu_cai_analysis_plots\' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v\configurations_graphs'];


%File path when working on MaCOS
path = ['hu_cai_analysis_plots/' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v/configurations_graphs'];

mkdir(path)

% [H,~,~,~,~,~,~,~,~,~] = hu_cai_model_with_fixed_sinks_and_sources(size_graph_x, size_graph_y,source_loc,sinks,source_bound,tspan,c_0,weights);
% 
% for i = 1:size(ordered_unique_k_configs,1)
%                 figure
%                 HWidths = 10*ordered_unique_k_configs(i,:)/max(ordered_unique_k_configs(i,:));
%                 num_nodes = numnodes(H);
%                 p = plot(H,'Layout','force','LineWidth',HWidths);
%                 p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
%                 saveas(gca, fullfile(path4,['config' int2str(i) '.fig']), 'fig');
%                 saveas(gca, fullfile(path4,['config' int2str(i) '.jpeg']), 'jpeg')
% end


%File path when working on Windows
%path2 = ['hu_cai_raw_data_plots\' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v'];

%File path when working on MacOS
path2 = ['hu_cai_raw_data_plots/' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v'];

mkdir(path2);


figure('Name','Number of Configurations vs. Edges Used')
plot(absolute_A,configs_num);
title('Number of Configurations vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Number of Configurations')
saveas(gca, fullfile(path2,'num_configs_v_number_of_edges_used.fig'), 'fig');
saveas(gca, fullfile(path2,'num_configs_v_number_of_edges.jpeg'), 'jpeg')

figure('Name','Total Frequency vs Edges Used')
plot(absolute_A, total_frequency)
title('Total Frequency vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Total Frequency')
saveas(gca, fullfile(path2,'freq_of_configs_v_number_of_edges_used.fig'), 'fig');
saveas(gca, fullfile(path2,'freq_of_configs_v_number_of_edges_used.jpeg'), 'jpeg')



top_freqs = 5;
[max_freqs, max_freqs_idx] = maxk(configuration_frequency,top_freqs);
% [max_dists, max_dists_idx] = maxk(avg_distance_array,top_freqs);
[~,max_freq_idx] = max(configuration_frequency);
%[~,min_freq_idx] = min(configuration_frequency);
% avg_value = mean(configuration_frequency);
% [common_value, common_value_idx] = min(abs(configuration_frequency-avg_value));
% tri_configs = [ordered_unique_k_configs(max_freq_idx,:); ordered_unique_k_configs(max_freqs_idx(3),:); ordered_unique_k_configs(max_freqs_idx(4),:)];

top_unique_k_configurations = ordered_unique_k_configs(max_freqs_idx,:);
top_configurations = ordered_configs(max_freqs_idx,:);



seen_configs = zeros(top_freqs,size(total_configs,2));
%File path when working on Windows
% path3 = ['hu_cai_analysis_plots\' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'mean_' num2str(v) 'variance\top_configs2by2'];

%File path when working on MacOS
path3 = ['hu_cai_analysis_plots/' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v/top_configs2by2'];

mkdir(path3)
for i = 1:size(top_configurations,1)
    h = 0;
    while h == 0
        [H,~,~,~,Q,~,k,~,~,~] = hu_cai_model_with_fixed_sinks_and_sources(size_graph_x, size_graph_y,source_loc,sinks,source_bound,tspan,c_0,weights,m,v);
        c = (abs(Q')>0);
        if c == top_configurations(i,:)
            if ismember(c,seen_configs,'rows') == 0
                h = h + 1;
                seen_configs(i,:) = (abs(Q)>0)';
                figure
                HWidths = 10*k/max(k);
                num_nodes = numnodes(H);
                p = plot(H,'Layout','force','LineWidth',HWidths);
                p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);
                saveas(gca, fullfile(path3,['top_config' int2str(i) '.fig']), 'fig');
                saveas(gca, fullfile(path3,['top_config' int2str(i) '.jpeg']), 'jpeg')
                
            end
        end
    end
end

%Windows file path
% filename = ['data\analysis_hu_cai' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v.mat'];

%MacOS file path
filename = ['data/analysis_hu_cai' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'm_' num2str(v) 'v.mat'];


% save(filename, 'distance_matrix', 'ordered_configs_idx', 'ordered_configs', 'total_frequency', 'configuration_frequency', 'top_configurations', 'tri_configs', 'top_unique_k_configurations', 'ordered_unique_k_configs')
save(filename,'ordered_configs_idx', 'ordered_configs', 'total_frequency', 'configuration_frequency', 'top_configurations', 'top_unique_k_configurations', 'ordered_unique_k_configs')
%%

%Windows file path
% path4 = ['hu_cai_analysis_plots\' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'mean_' num2str(v) 'variance\by_geometry'];

%MacOS file path
path4 = ['hu_cai_analysis_plots/' int2str(size_graph_x) '_by_' int2str(size_graph_y) '_' int2str(length(source_loc)) 'source_' int2str(length(sinks)) 'sinks_' num2str(m) 'mean_' num2str(v) 'variance/by_geometry'];

mkdir(path4)
f2 = figure('Name','Energy Dissapation vs. Number of Edges Used');
scatter(total_used_edges, total_energy_array)
title('Energy Dissapation vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Energy Dissapation')
saveas(gca, fullfile(path4,'energy_edges.fig'), 'fig');
saveas(gca, fullfile(path4,'energy_edges.jpeg'), 'jpeg')

figure('Name','Configuration Frequency');
plot(ordered_configs_idx, configuration_frequency)
title('Total Frequency vs. Configuration')
xlabel('Configuration')
ylabel('Total Frequency')
for i = 1:length(configs_num)
    vline(sum(configs_num(1:i)),'r')
end
    
saveas(gca, fullfile(path4,'total_cofiguration_freq_v_configuration.fig'), 'fig');
saveas(gca, fullfile(path4,'total_cofiguration_freq_v_configuration.jpeg'), 'jpeg')

% figure('Name','Average Distance From Other Optimum');
% plot(ordered_configs_idx, avg_distance_array)
% title('Average Distance From Other Optimum vs. Configuration')
% xlabel('Configuration')
% ylabel('Average Distance From Other Optimum')
% saveas(gca, fullfile(path1,'avg_dist.fig'), 'fig');
% saveas(gca, fullfile(path1,'avg_dist.jpeg'), 'jpeg');


figure('Name','Energy at Each Configuration');
scatter(energy_array_idx, total_energy_array)
title('Energy at Each Configuration vs. Configuration')
xlabel('Configuration')
ylabel('Energy at Each Configuration')

hold on
plot(ordered_configs_idx, avg_energy_array)

hold off
saveas(gca, fullfile(path4,'config_energy'), 'fig');
saveas(gca, fullfile(path4,'config_energy'), 'jpeg');


% figure('Name', 'Avg Energy Distance at Each Configuration');
% plot(ordered_configs_idx, avg_energy_difference)
% title('Avg Energy Distance at Each Configuration vs. Configuration')
% xlabel('Configuration')
% ylabel('Avg Energy Distance at Each Configuration')
% saveas(gca, fullfile(path1,'avg_energy_dist'), 'fig');
% saveas(gca, fullfile(path1,'avg_energy'), 'jpeg');

figure('Name', 'Energy vs Configuration Frequency');
scatter(big_E, log(configuration_frequency))
%set(gca, 'yscale', 'log')
title('Energy at Configuration vs Configuration Frequency')
xlabel('Configuration Energy')
ylabel('Configuration Frequency')
saveas(gca, fullfile(path4,'energy_v_config_frequency'), 'fig');
saveas(gca, fullfile(path4,'energy_v_config_frequency'), 'jpeg');

% figure('Name', 'Energy vs Sum Configuration Frequency');
% scatter(little_E, avg_freq_for_energy_level)
% title('Energy at Configuration vs Configuration Frequency')
% xlabel('Configuration Energy')
% ylabel('Configuration Frequency')
% saveas(gca, fullfile(path4,'energy_v_sum_frequency'), 'fig');
% saveas(gca, fullfile(path4,'energy_v_sum_frequency'), 'jpeg');

figure('Name', '(Binned)Energy vs Configuration Frequency');
scatter(avg_of_binned_energy, log(avg_of_binned_configuration_freq))
%set(gca, 'yscale', 'log')
title('Energy at Configuration vs Configuration Frequency')
xlabel('Configuration Energy')
ylabel('Configuration Frequency')
saveas(gca, fullfile(path4,'binned_energy_v_config_frequency'), 'fig');
saveas(gca, fullfile(path4,'binned_energy_v_config_frequency'), 'jpeg');


% figure('Name', 'Histogram of Frequency of Frequencies');
% x = log10(configuration_frequency);
%  [~,edges] = histcounts(x);
%  y = 10.^edges;
% histogram(configuration_frequency,y)
% set(gca, 'xscale','log')
% title('Number of Different Configurations vs Frequency of Configuration Frequency')
% xlabel('Frequency of Configuration Frequency')
% ylabel('Number of Different Configurations')
% saveas(gca, fullfile(path1,'freq_histogram'), 'fig');
% saveas(gca, fullfile(path1,'freq_histogram'), 'jpeg');

% figure('Name', 'Number of Edges Used vs Configuration');
% plot(ordered_configs_idx,edges_array)
% saveas(gca, fullfile(path1,'freq_histogram'), 'fig');
% saveas(gca, fullfile(path1,'freq_histogram'), 'jpeg');


% figure('Name', 'Minimum Distance from Optimum vs Configuration');
% scatter(ordered_configs_idx, min_distance)
% title('Minimum Distance from Optimum vs Configuration')
% xlabel('Configuration')
% ylabel('Distance')
% saveas(gca, fullfile(path1,'min_dist'), 'fig');
% saveas(gca, fullfile(path1,'min_dist'), 'jpeg');
% 
% figure('Name', 'Closest Configuration vs Configuration');
% scatter(config_idx_min, min_idx)
% title('Closest Configuration vs Configuration')
% xlabel('Configuration')
% ylabel('Closest Configuration')
% saveas(gca, fullfile(path1,'closest_configuration'), 'fig');
% saveas(gca, fullfile(path1,'closest_configuration'), 'jpeg');
% 
% figure('Name', 'Frequency vs. Closest Configurations');
% scatter(ordered_configs_idx, min_freq)
% title('Frequency vs. Closest Configurations')
% xlabel('Closest Configuration')
% ylabel('Frquency')
% saveas(gca, fullfile(path1,'closest_freq'), 'fig');
% saveas(gca, fullfile(path1,'closest_freq'), 'jpeg');
% 
% 
% figure('Name', 'Histogram of Frequency of Minimum Distance Frequencies');
% x = log10(min_freq);
%  [~,edges] = histcounts(x);
%  y = 10.^edges;
% histogram(min_freq,y)
% set(gca, 'xscale','log')
% title('Number of Different Configurations vs. Frequency of Minimum Distance Frequencies')
% xlabel('Frequency of Minimum Distance Frequencies')
% ylabel('Number of Different Configurations')
% saveas(gca, fullfile(path1,'closest_histogram'), 'fig');
% saveas(gca, fullfile(path1,'closest_histogram'), 'jpeg');
% 
% figure('Name', 'Maximum Distance from Optimum vs Configuration');
% scatter(ordered_configs_idx, max_distance)
% title('Maximum Distance from Optimum vs Configuration')
% xlabel('Configuration')
% ylabel('Distance')
% saveas(gca, fullfile(path1,'max_dist'), 'fig');
% saveas(gca, fullfile(path1,'max_dist'), 'jpeg');
% 
% 
% figure('Name', 'Farthest Configuration vs Configuration');
% scatter(config_idx_max, max_idx)
% title('Farthest Configuration vs Configuration')
% xlabel('Configuration')
% ylabel('Farthest Configuration')
% saveas(gca, fullfile(path1,'farthest_configuration'), 'fig');
% saveas(gca, fullfile(path1,'farthest_configuration'), 'jpeg');
% 
% 
% figure('Name', 'Frequency vs. Farthest Configurations');
% scatter(ordered_configs_idx, max_freq)
% title('Frequency vs. Farthest Configurations')
% xlabel('Farthest Configuration')
% ylabel('Frquency')
% saveas(gca, fullfile(path1,'farthest_freq'), 'fig');
% saveas(gca, fullfile(path1,'farthest_freq'), 'jpeg');
% 
% 
% figure('Name', 'Histogram of Frequency of Maximum Distance Frequencies');
% x = log10(max_freq);
%  [~,edges] = histcounts(x);
%  y = 10.^edges;
% histogram(max_freq,y)
% set(gca, 'xscale','log')
% title('Number of Different Configurations vs. Frequency of Maximum Distance Frequencies')
% xlabel('Frequency of Maximum Distance Frequencies')
% ylabel('Number of Different Configurations')
% saveas(gca, fullfile(path1,'farthest_histogram'), 'fig');
% saveas(gca, fullfile(path1,'farthest_histogram'), 'jpeg');



end


