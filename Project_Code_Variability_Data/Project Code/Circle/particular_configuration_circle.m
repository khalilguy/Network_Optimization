function [configuration_frequency,ordered_configs,fg,fg2,fg3,fg4,fg5,fg6] = particular_configuration_circle(source_loc,sinks,source_bound,radius,trials,target)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tol = .001;
load possible_configurationsr4.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array trials_array diff_array config_array

total_edgelist = (abs(total_Q_matrix) > 0)';
ordered_configs = zeros(size(total_configs,1), size(total_configs,2));
configuration_frequency = zeros(size(total_configs,1),1);
energy_array_idx = zeros(total_trials,1);
avg_energy_array = zeros(size(total_configs,1), 1);
m = 0;
for j = 1:size(absolute_A ,1)
    for i = 1:size(total_configs,1)
        if absolute_A(j,1) == sum(total_configs(i,:))
            m = m + 1;
            ordered_configs(m,:) = total_configs(i,:);
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

ordered_configs_idx = 1:size(ordered_configs,1);
edges_array = zeros(current_config,1);
distance_idx_array = zeros((current_config - 1)*current_config, 1);
distance_array = zeros((current_config - 1)*current_config, 1);
avg_distance_array  = zeros(current_config,1);
avg_energy_difference = zeros(current_config,1);
k = 1;



for i = 1:current_config 
    p = 0;
    config_i = ordered_configs(i,:);
    
    for g = 1:current_config
        if i~=g
            avg_energy_difference(i,1) = avg_energy_difference(i,1) + sqrt((avg_energy_array(i,1) - avg_energy_array(g,1))^2);
        end
    end
    
    while p < current_config
        
        p = p + 1;
        if p ~= i
            config_j = ordered_configs(p,:);
            different = 0;
            for j = 1:size(config_i,2)
                if config_i(1,j) ~= config_j(1,j)
                    avg_distance_array(i,1) = avg_distance_array(i,1) + 1;
                    different = different + 1;
                end
            end
            distance_idx_array(k,1) = i;
            distance_array(k,1) = different;
            
            k = k + 1;
        end
    end
    avg_distance_array(i,1) = avg_distance_array(i,1)/(current_config - 1);
    avg_energy_difference(i,1) = avg_energy_difference(i,1)/(current_config-1);
    edges_array(i,1) = sum(config_i);
end



top_freqs = 4;
[max_freqs, max_freqs_idx] = maxk(configuration_frequency,top_freqs);
[max_dists, max_dists_idx] = maxk(avg_distance_array,top_freqs);
top_configurations = zeros(top_freqs,size(total_configs,2));


for i = 1:top_freqs
    top_configurations(i,:) = ordered_configs(max_freqs_idx(i,1),:);
end

h = 0;

seen_configs = zeros(top_freqs,size(total_configs,2));

for i = 1:top_freqs
    h = 0;
    while h == 0
    [Q,k,H] = model_on_circle(source_loc,sinks,source_bound,radius,tol);
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
            end
        end
    end
end
    
set(0,'DefaultFigureVisible','on');
fg = figure('Name','Configuration Frequency');
bar(ordered_configs_idx, configuration_frequency)
title('Total Frequency vs. Configuration')
xlabel('Configuration')
ylabel('Total Frequency')

fg2 = figure('Name','Average Distance From Other Optimum');
bar(ordered_configs_idx, avg_distance_array)
title('Average Distance From Other Optimum vs. Configuration')
xlabel('Configuration')
ylabel('Average Distance From Other Optimum')

fg3 = figure('Name','Energy at Each Configuration');
scatter(energy_array_idx, total_energy_array)
title('Energy at Each Configuration vs. Configuration')
xlabel('Configuration')
ylabel('Energy at Each Configuration')

hold on
plot(ordered_configs_idx, avg_energy_array)

hold off

fg4 = figure('Name', 'Avg Energy Distance at Each Configuration');
bar(ordered_configs_idx, avg_energy_difference)
title('Avg Energy Distance at Each Configuration vs. Configuration')
xlabel('Configuration')
ylabel('Avg Energy Distance at Each Configuration')

fg5 = figure('Name', 'Histogram of Frequency of Frequencies');
histogram(configuration_frequency,15)
title('Frequency of Configuration Frequency vs. Number of Different Configurations')
xlabel('Frequency of Configuration Frequency')
ylabel('Number of Different Configurations')

fg6 = figure('Name', 'Histogram of Frequency of Frequencies');
plot(ordered_configs_idx,edges_array)
end

