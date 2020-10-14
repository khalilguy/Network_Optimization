function [H,time,C,k_i,Q_final,Q_initial,k_final, Q_flow_direction_final,changed_flow_direction]= hu_cai(tspan,c_0,s,t,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions,weights,m,v)

num_edges = length(s);
weights = lognrnd(m,v,1,num_edges);
c_0  = sum(weights.^1/2);
[k_i,I,it,F,p,num_nodes,H] = create_digraph(s,t,weights,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions);
Q_initial = calculate_flows(k_i,I,it,F,p,num_nodes);
positive  = Q_initial > 0;
negative = (Q_initial < 0)*-1;
Q_flow_direction_initial = positive + negative;
[time,C] = ode15s(@(t,c)adaptation_ode(t,c,c_0,I,it,F,p,num_nodes),tspan, k_i);
k_final = C(end,:)';
k_final = max(k_final,10^-6);
Q_final = calculate_flows(k_final,I,it,F,p,num_nodes);

positive  = Q_final > 0;
negative = (Q_final < 0)*-1;
Q_flow_direction_final = positive + negative;

changed_flow_direction = 0;

for i = 1:length(Q_final)
    if Q_flow_direction_initial(i,1) ~= Q_flow_direction_final(i,1)
        changed_flow_direction = 1;
        break
    end
end

for i = 1:size(Q_final)
    if abs(Q_final(i,1)) < sum(source_boundary_conditions)/length(Q_final)
        Q_final(i,1) = 0;
    end
end

H.Edges.Weight = k_final;


end


% num_edges = length(s);
% source_loc = 1;
% sink_loc = [8,13]';
% source_boundary_conditions = 500;
% sink_boundary_conditions = [-250,-250]';
% weights = ones(1,num_edges) + .5*rand(1,num_edges);
% [k,I,it,F,p,num_nodes,H] = create_digraph(s,t,weights,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions);
% Q = calculate_flows(k,I,it,F,p,num_nodes);
% HWidths = 5*k/max(k);
% c_0 = 1;
% [t,C] = ode45(@(t,c)adaptation_ode(t,c,c_0,I,it,F,p,num_nodes),tspan, k);
% 

% figure(2)
% plot(t, C(:,:))
% figure(3)
% HWidths = 5*C(end,:)./max(C(end,:));
% Q = calculate_flows(C(end,:),I,it,F,p,num_nodes);
% plot(H,'Layout','force','EdgeLabel',Q,'LineWidth',HWidths)

function dCdt = adaptation_ode(t,c,c_0,I,it,F,p,num_nodes)

    gamma = 1/2;
    tau = gamma*c_0;
    c = max(c,10^-6);
    q = calculate_flows(c,I,it,F,p,num_nodes);
    dCdt = ((q.^2)./(c.^(gamma + 1)) - ones(length(q),1)*tau).*c;
    dCdt(c<10^-6) = max(dCdt(c<10^-6),0);
end

