%%
%0,1
source_loc = [8,4]';
sinks = [7,1]';

% 2,2
% source_loc = [32 36 40]';
% sinks = [1,7,11]';

%2,3
% source_loc = [15,11,7,1]';
% sinks = [52,48,44,40]';


%2,9
% source_loc = [27,23,19,15]';
% sinks = [112,108,104,100]';

%R = 1
% source_loc = 2
% sinks = [1,3,4]';
%R = 4
source_loc = 17;
sinks = [27,40,5]';

target = 150;
radius = 4;
trials = 1000;
source_bound = 500;

[A,f1,f2] = energy_and_frequency_on_circle(source_loc,sinks,source_bound,trials,radius);
% [configuration_frequency,ordered_configs,fg,fg2,fg3,fg4,fg5,fg6] = particular_configuration(2,2,source_loc,sinks,source_bound);
% [A2,f3,f4] = energy_and_frequency_vs_edges(0,1,source_loc,sinks,source_bound,trials);
% [total_configs,config_nums,total_trials,absolute_A,total_frequency] = possible_configurations(0,1,source_loc,sinks,source_bound,trials,target);












