function [current_config,configs_num,total_configs,total_trials,absolute_A,total_frequency,total_Q_matrix,total_energy_array] = possible_configurations(size_graph_x, size_graph_y,source_loc,sinks,source_bound,trials,target)

%parameters
%possible_configuraions - Function used to find all the possible local
%optima for a particular layout of a network. Also finds the local optima
%and maintains the direction of flow in the data

%size_graph_x - provides width of graph

%size_graph_y - provides length of graph

%source_loc - location(s) of souce(s)

%sinks - location(s) of sink(s)
%source_bound - the source boundary conditions for the network. This
%boundary condition is split evenly among the number of sources and the
%number of sinks to respect laws of flow conservation
%trials - the search for the local optima batches searches based on the
%number of trials. Higher values for trials may imporve the number of
%networks it able to find but it will take more time

%target - as the algorithm finds more local optimum, the search will begin
%to asymptote towards some value. If the difference between the previous number 
%of trials and the current number of trials is zero, then the search has
%begun to become asymptotic. But since the search is random there
%could still be some values it has not found. Target tells how many times
%you want to find a difference of zero. Higher values mean you will find
%more local optimum, lower means you may still have some left to find

%returns, as well as saves in a Matlab data file
%current_config - the number of local optima the algorithm found
%configs_num - an array where each element represents the number of local
%optima found for a particular number of edges. It is ordered in the same
%ways as absolute_A

%total_configs - a logical matrix consisting of local optima. Each row
%represents a local optima and the colums represent an edge along the
%lattice. 1 represents non-zero flow across and edge, 0 represents zero
%flow across an edge

%total_trials - the total number of trials run

%absolute_A - an array where each element represents the possible number of
%edges that could be used. This array is sorted from lowest to highest

%total_frequency - an array where each element represents the frequency at
%which a local optima which used a certain number of edges. For the (i,1)
%element in total_frequency, the number of edges for that element is the
%(i,1) element in absolute_A

%total_Q_matrix - like total_configs but gives all the final configurations
%over all the trials run. Also, it fully represents the local optimum where each 
%row represents an edge and the element represent the flow along that edges.
%The number of columns corresponds to the total_trials variable

%total_energy_array - the energy dissipation correpondint the (:,i)
%local_optima in total_Q_matrix

%This function can take a while to run depending on the size of the network
%Another function uses alot of the variables from this, so they are saved 
%to a Matlab data file to reduce the number of times they need to be recalculated

%Other saved variable:
%Any variable with the word "flow" in it is used represent the above
%variables while still maintaining the direction of flow in the data.
%Depending on how you choose to inititally represent the array for this,
%you may get different values for the number of trials run
%total_flow_trials  - total number of trials run for the finding the
%possible configurations for flow direction
%changed_direction_array - a logical array where each element tells if
%after the first iteration of the relaxation algorithm, the flow across adn
%edge changed direction
%total_flow_configs -total number of local optima found where a local
%optima is represented as by the flows and its direction 
%trials_array - each element represents the number of trials that have been
%run. They are spaced by the variable trials

%config_array - an array that shows the number of local optima that have
%been found after so many trials have been run. If the (i,1) element in
%config_array shows the number of local optima that have been found at a
%number of trials, trials_array(i,1) is the number of trials
%flow_trials_array - same thing as trials array, but for the flow
%calculations

[Q,~,~] = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
total_configs = logical.empty(0,size(Q,1));
absolute_A = double.empty(0,1);
current_config = 100;
prev_config = 0;
prev_flow_config = 0;
numit = 1;
numit_flow = 1;
total_trials = 0;
total_flow_trials =0;
found = 0;
found_flow = 0;
gamma = .5;

while found < target && found_flow < target
    for i = 1:trials
        [Q,~,~,changed_direction,~,~,Qf] = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
        used_edges = 0;
        E = 0;
        %calculating energy dissipation
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

        %saving all the important information each trial
        total_Q_matrix(:,i+total_trials) = Q;
        total_used_edges(i + total_trials,1) = used_edges;
        total_energy_array(i+total_trials,1) = E;
        total_Qf_matrix(:,i+total_flow_trials) = Qf;
        changed_direction_array(i+total_trials,1) = changed_direction;
    end
    
    %puts the local optima purely in terms of whether or not they have a 
    total_edgelist = abs(total_Q_matrix)> 0;
    
    [total_configs,idxed,idxtc] = unique(total_edgelist', 'rows');
    total_flow_configs = total_Q_matrix(:,idxed)';
    
    
    current_config = size(total_configs,1)
    current_config_flow = size(total_flow_configs,1)
    
    diff = sqrt((current_config - prev_config)^2)
    diff_flow = sqrt((current_config_flow - prev_flow_config)^2)
    
    prev_config = current_config;
    prev_flow_config = current_config_flow;
    
    total_trials = total_trials + trials;
    total_flow_trials = total_flow_trials + trials;
    
    flows_config_array(numit_flow,1) = current_config_flow;
    flow_trials_array(numit_flow,1) = total_flow_trials;
    config_array(numit,1) = current_config;
    trials_array(numit,1) = total_trials;
    diff_array(numit,1) = diff;
    
    numit = numit + 1
    numit_flow = numit_flow + 1
  
    if diff == 0
        found = found + 1
    end
    
    if diff_flow == 0
        found_flow = found_flow + 1
    end

    absolute_A = unique(total_used_edges);
    
end


while found < target
    gamma = .5;
    for i = 1:trials
        Q = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
        used_edges = 0;
        E = 0;
        for j = 1: size(Q)
            E = E + (Q(j,1)^2)^(gamma/(gamma +1));
            if abs(Q(j,1)) > 0
                used_edges = used_edges + 1;
            end
        end

        total_Q_matrix(:,i+total_trials) = Q;
        total_used_edges(i + total_trials,1) = used_edges;
        total_energy_array(i+total_trials,1) = E;
    end

    total_edgelist = abs(total_Q_matrix)> 0;

    total_configs = unique(total_edgelist', 'rows');

    current_config = size(total_configs,1)
    diff = sqrt((current_config - prev_config)^2)
    prev_config = current_config;
    total_trials = total_trials + trials;
    config_array(numit,1) = current_config;
    trials_array(numit,1) = total_trials;
    diff_array(numit,1) = diff;
    numit = numit + 1

    absolute_A = unique(total_used_edges);

    if diff == 0
        found = found + 1
    end
end



while found_flow < target
    for i = 1:trials
        [~,~,~,changed_direction,~,~,Qf] = model_with_fixed_sinks_and_sources(size_graph_x,size_graph_y,source_loc,sinks,source_bound);
        total_Qf_matrix(:,i+total_flow_trials) = Qf;
        changed_direction_array(i+total_flow_trials,1) = changed_direction;
    end


    total_flow_configs = unique(total_Qf_matrix','rows');
    current_config_flow = size(total_flow_configs,1)
    diff_flow = sqrt((current_config_flow - prev_flow_config)^2)
    prev_flow_config = current_config_flow;
    total_flow_trials = total_flow_trials + trials;
    numit_flow = numit_flow + 1


    if diff_flow == 0
        found_flow = found_flow + 1
    end
end

    

absolute_A = sort(absolute_A);
configs_num = zeros(length(absolute_A),1);

for i = 1:length(absolute_A)
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


path = 'possible_plots';
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
saveas(gca, fullfile(path,'trials_diff_configs.fig'), 'fig');
saveas(gca, fullfile(path,'trials_diff_configs.jpeg'), 'jpeg');

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
saveas(gca, fullfile(path,'num_freq_configs.fig'), 'fig');
saveas(gca, fullfile(path,'num_freq_configs.jpeg'), 'jpeg')

% save possible_configurations0by1.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges changed_direction_array total_flow_configs
save possible_configurations2by2.mat current_config current_config_flow configs_num total_configs total_trials total_flow_trials absolute_A total_frequency total_Q_matrix total_Qf_matrix total_energy_array total_used_edges changed_direction_array total_flow_configs trials_array config_array flow_trials_array
% save possible_configurations2by3.mat current_config configs_num total_configs total_trials absolute_A total_frequency total_Q_matrix total_energy_array total_used_edges changed_direction_array total_flow_configs




end

