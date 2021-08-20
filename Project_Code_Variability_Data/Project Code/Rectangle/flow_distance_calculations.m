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