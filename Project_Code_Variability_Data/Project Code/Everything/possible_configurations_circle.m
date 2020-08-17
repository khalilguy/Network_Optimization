function [current_config,configs_num,total_configs,total_trials,absolute_A,total_frequency,total_Q_matrix,total_energy_array] = possible_configurations_circle(source_loc,sinks,source_bound,radius,tol,trials,target)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Q,~,~] = model_on_circle(source_loc,sinks,source_bound,radius,tol);
total_configs = logical.empty(0,size(Q,1));
absolute_A = double.empty(0,1);
current_config = 100;
prev_config = 0;
num_it = 1;
total_trials = 0;
found = 0;

while found < target
    used_edges_array = zeros(trials,1);
    energy_array = zeros(trials,1);
    gamma = .5;
    [Q,~,~] = model_on_circle(source_loc,sinks,source_bound,radius,tol);
    Q_matrix = zeros(size(Q,1),trials);
    for i = 1:trials
            [Q,~,~] = model_on_circle(source_loc,sinks,source_bound,radius,tol);
            used_edges = 0;
            E = 0;
            for j = 1: size(Q)
                E = E + (Q(j,1)^2)^(gamma/(gamma +1));
                if abs(Q(j,1)) > 0
                    used_edges = used_edges + 1;
                end
            end
    %         if used_edges == 25
    %             disp('found')
    %             break
    %         end
            Q_matrix(:,i) = Q;
            total_Q_matrix(:,i+total_trials) = Q;
            used_edges_array(i,1) = used_edges;
            total_used_edges(i + total_trials,1) = used_edges;
            total_energy_array(i+total_trials,1) = E;
            energy_array(i,1) = E;
    end

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

    edgelist = abs(Q_matrix) > 0;

    new_freq = frequency_array;
    new_used_edges = used_edges_array;

    for i = 1:size(A)
        edges = A(i,1);
        config_matrix = zeros(size(Q,1),frequency_array(i,1));
        j = 1;
        
        for k = 1:trials
            if edges == used_edges_array(k,1)
                config_matrix(:,j) = edgelist(:,k);
                j = j + 1;
            end
        end

        config_matrix = unique(config_matrix', 'rows');
        if size(total_configs,1) == 0
            total_configs = config_matrix;
        elseif size(total_configs,1) > 0
            for h = 1:size(config_matrix,1)
                if ismember(config_matrix(h,:), total_configs, 'rows') == 0
                    total_configs = [total_configs; config_matrix(h,:)];
                end
            end
        end
    end
    
    if size(absolute_A,1) == 0
            absolute_A = A;
    elseif size(absolute_A,1) > 0
        for h = 1:size(A,1)
            if ismember(A(h,:), absolute_A, 'rows') == 0
                absolute_A = [absolute_A; A(h,:)];
            end
        end
    end
    current_config = size(total_configs,1)
    diff = sqrt((current_config - prev_config)^2)
    prev_config = current_config;
    
    total_trials = total_trials + trials
    config_array(num_it,1) = current_config;
    trials_array(num_it,1) = total_trials;
    diff_array(num_it,1) = diff;
    num_it = num_it + 1
    
    if diff == 0
        found = found + 1
    end
    
end


absolute_A = sort(absolute_A);
configs_num = zeros(size(absolute_A,1),1);

for i = 1:size(absolute_A)
    for j = 1:size(total_configs)
        if sum(total_configs(j,:)) == absolute_A(i,1)
            configs_num(i,1) = configs_num(i,1) + 1;
        end
    end  
end


total_frequency = zeros(size(absolute_A,1),1);
for i = 1:total_trials
        for j = 1:size(absolute_A,1)
            if total_used_edges(i,1) == absolute_A(j,1)
                total_frequency(j,1) = total_frequency(j,1) + 1;
            end
        end
end



set(0,'DefaultFigureVisible','on');
fg = figure('Name','Number of Trials Performed vs. Distance from Total Configurations');
subplot(1,2,1)
plot(trials_array, diff_array);
title('Number of Trials Performed vs. Distance from Total Configurations')
xlabel('Number of Trials Performed')
ylabel('Distance from Total Configurations')

subplot(1,2,2)
plot(trials_array,config_array);
title('Number of Trials Performed vs. Number of Configurations Found')
xlabel('Number of Trials Performed')
ylabel('Number of Configurations Found')

fg2 = figure('Name','Number of Configurations for a Particular Number of Edges Used vs. Number of Edges Used');
subplot(1,2,1)
plot(absolute_A,configs_num)
title('Number of Configurations for a Particular Number of Edges Used vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Number of Configurations')

subplot(1,2,2)
plot(absolute_A, total_frequency)
title('Total Frequency vs. Number of Edges Used')
xlabel('Number of Edges Used')
ylabel('Total Frequency')

save possible_configurationsr4.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array trials_array diff_array config_array
end

