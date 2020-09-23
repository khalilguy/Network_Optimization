tspan = 0:100;
c_0 = 100;
weights = 1;
target = 10;
trials = 1000;
source_bound = 500;
v = 1;
m = -v.^2/2;

source_loc = [8,4]';
sinks = [7,1]';

generating_data_hu_cai(0,1,source_loc,sinks,source_bound,trials,target,tspan,c_0,weights,m,v);
