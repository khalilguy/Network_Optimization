
%R = 1
% source_loc = 2
% sinks = [1,3,4]';
%R = 4
source_loc = 17;
sinks = [27,40,5]';

target = 10;
radius = 4;
trials = 1500;
source_bound = 500;

% [A,f1,f2] = energy_and_frequency_on_circle(source_loc,sinks,source_bound,trials,radius);
%  [total_configs,config_nums,total_trials,absolute_A,total_frequency] = possible_configurations_circle(source_loc,sinks,source_bound,radius,trials,target);
[configuration_frequency,ordered_configs,fg,fg2,fg3,fg4,fg5,fg6] = particular_configuration_circle(source_loc,sinks,source_bound,radius,trials,target);
