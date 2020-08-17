function [final_configs, initial_distance, I] = initial_k(size_graph_x, size_graph_y,source_loc,sinks,source_bound,trials)

%Objective of function: Find n number of initial conductances which lead to
%a particular number of edges used in the network
%   Detailed explanation goes here

[~,k,~,~,k_i,edgeidx] = model_with_fixed_sinks_and_sources(size_graph_x ,size_graph_y,source_loc,sinks,source_bound);
used_edges_array = zeros(trials,1);
initialK = zeros(length(k),trials);

load particular_configuration2by2.mat distance_matrix ordered_configs
load particular_flow_configuration2by2.mat total_flow_configs ordered_flow_configs

initial_distance = zeros(trials,size(ordered_flow_configs,1));
final_configs = zeros(trials,1);

for i = 1:trials
    [Q,k,~,~,k_i,~,~,Q_i] = model_with_fixed_sinks_and_sources(size_graph_x ,size_graph_y,source_loc,sinks,source_bound);
    for j = 1:size(ordered_flow_configs)
         initial_distance(i,j) = sqrt(sum((Q_i - ordered_flow_configs(j,:)').^2));%Euclidean distance may not be the best metric
         if (abs(Q)>0)' == ordered_configs(j,:)
             final_configs(i,1) = j;    
         end
    end
end

[~,I] = min(initial_distance,[],2);


figure
plot(1:trials, I)
figure
plot(1:trials,final_configs)
 
freq_closest_config = zeros(size(ordered_flow_configs,1),1);
freq_final_config = zeros(size(ordered_flow_configs,1),1);

for i = 1: size(ordered_flow_configs,1)
    for j = 1:trials
        if i == I(j,1)
            freq_closest_config(i,1) = freq_closest_config(i,1) + 1;
        end
        if i == final_configs(j,1)
            freq_final_config(i,1) = freq_final_config(i,1) + 1;
        end
    end
end

figure
plot(1:length(freq_closest_config),freq_closest_config )
figure
plot(1:length(freq_closest_config),freq_final_config)            


%  for i = 1:trials
%     [Q,k,~,~,k_i,~,Q_i] = model_with_fixed_sinks_and_sources(size_graph_x ,size_graph_y,source_loc,sinks,source_bound);
%     used_edges = 0;
%     for j = 1: size(Q)
%         if abs(Q(j,1)) > 0
%             used_edges = used_edges + 1;
%         end
%     end
%     used_edges_array(i,1) = used_edges;
%     initialK(:,i)= k_i;
%  end
%  
%  
%  A = unique(used_edges_array);
%  frequency_array = zeros(length(A),1);
%  
%  for i = 1:length(A)
%      idx = 0;
%      for j = 1:trials
%          if A(i,1) == used_edges_array(j,1)
%              if idx == 0
%                  idx = j;
%              else
%                  idx = [idx j];
%              end
%              frequency_array(i,1) = frequency_array(i,1) + 1;
%          end
%      end
%      
%      k_it = zeros(length(k)*length(idx),1);
%      edges_idx = zeros(length(edgeidx)*length(idx),1);
%      
%      h = 1;
%      
%      for j = 1:length(k):length(k_it)
%          k_it(j:j+length(k)-1,1) = initialK(:,idx(1,h));
%          edges_idx(j:j+length(k)-1,1) = edgeidx;
%          h = h +1;
%      end
%  end
%                
         
end

