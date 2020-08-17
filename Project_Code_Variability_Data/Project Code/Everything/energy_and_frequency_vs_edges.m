function [A,f1,f2] = energy_and_frequency_vs_edges(size_graph_x, size_graph_y,source_loc,sinks,source_bound,trials)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

used_edges_array = zeros(trials,1);
energy_array = zeros(trials,1);
gamma = .5;

[Q,k,H] = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
Q_matrix = zeros(size(Q,1),trials);

for i = 1:trials
        [Q,k,H] = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
        used_edges = 0;
        E = 0;
        for j = 1: size(Q)
            E = E + (Q(j,1)^2)^(gamma/(gamma +1));
            if abs(Q(j,1)) > 0
                used_edges = used_edges + 1;
            end
        end
        if used_edges == 31
            disp('found')
            break
        end
        Q_matrix(:,i) = Q;
        used_edges_array(i,1) = used_edges;
        energy_array(i,1) = E;
end



f1 = figure('Name', 'Graph');
HWidths = 10*k/max(k);
num_nodes = numnodes(H);
 p = plot(H,'Layout','force','LineWidth',HWidths);
p.NodeLabel = arrayfun(@num2str, 1:num_nodes, 'UniformOutput', false);


f2 = figure('Name','Energy Dissapation vs. Number of Edges Used');
subplot(1,3,1)
scatter(used_edges_array, energy_array)
title('Energy Dissapation vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Energy Dissapation')
A = unique(used_edges_array);
frequency_array = zeros(size(A,1),1);

for i = 1:trials
    edges = used_edges_array(i,1);
    for j = 1:size(A,1)
        if edges == A(j,1)
            frequency_array(j,1) = frequency_array(j,1) + 1;
        end
        
    end
end

subplot(1,3,2)
plot(A,frequency_array)
title('Frequency vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Frequency')

size(A,1)
edgelist = abs(Q_matrix) > 0;

configs_num = zeros(size(A,1),1);

new_freq = frequency_array;
new_used_edges = used_edges_array;
total_configs = logical.empty(0,size(Q,1));

for i = 1:size(A)
    edges = A(i,1);
    config_matrix = zeros(size(Q,1),frequency_array(i,1));
    j = 1;
    for k = 1:trials
        
        if edges == used_edges_array(k,1)
            
            config_matrix(:,j) = edgelist(:,k);
            j = j + 1;
%             used_edges_array(k,:) = [];
%             edgelist(:,k) = [];
%             trials = trials - 1;
        end
    end
    
    config_matrix = unique(config_matrix', 'rows');
    
    if size(total_configs,1) == 0
        total_configs = config_matrix;
    elseif size(total_configs,1) > 0
        for h = 1:size(config_matrix,1)
            con = config_matrix(h,:);
                if ismember(con, total_configs, 'rows') == 0
                    total_configs = [total_configs; con];
                end
        end
    end
    configs_num(i,1) = size(config_matrix,1);
end

size(total_configs)
frequency_array
configs_num
subplot(1,3,3)
plot(A,configs_num)

end

