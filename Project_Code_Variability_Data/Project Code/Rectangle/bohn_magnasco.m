function [Q,k,H,changed_direction,k_i,Qf,Q_i] = bohn_magnasco(s,t,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
num_edges = size(s,1);
weights = ones(1,num_edges) + .5*rand(1,num_edges);
[k,I,it,F,p,num_nodes,H] = create_digraph(s,t,weights,source_loc,sink_loc,source_boundary_conditions,sink_boundary_conditions);
Q = calculate_flows(k,I,it,F,p,num_nodes);

Q_i = Q;
k_i = k;
positive  = Q > 0;
negative = (Q < 0)*-1;
Q1 = positive + negative;

oldK = k;
my_step_size = .00001;
step_size = 100;
numit = 0;

sumK = 0;
%Calculate the conductance constraint
gamma = 1/2;
    for i = 1:length(k)
        sumK = sumK + k(i,1)^gamma;
    end

    bigK = sumK^(1/gamma);

while step_size > my_step_size
    %Reassigning Conductances
    sumQ = 0;

    %Calculating lambda
    for i = 1:length(Q)
        sumQ = sumQ + (Q(i,1)^2)^(gamma/(gamma + 1));
    end
    sumQ = sumQ^((gamma+1)/gamma);
    lambda = sumQ/(gamma*(bigK^(gamma + 1)));
    

    %Calculate the new conductance
     for i = 1:length(Q)
         if oldK(i,1) > 10^-5
            k(i,1) = ((Q(i,1)^2)/(lambda*gamma))^(1/(gamma +1));
         end
     end

     step_size = sqrt(sum((oldK - k).^2));
     oldK = k;


    %Recalculating Flows
    Q = calculate_flows(k,I,it,F,p,num_nodes);
    
    if numit == 0
        positive  = Q > 0;
        negative = (Q < 0)*-1;
        Q2 = positive + negative;
    
        changed_direction = 0;
        if Q1 ~= Q2
         changed_direction = 1;
        end
    end
    
    numit = numit + 1;

end

Qf = Q;
positive  = Qf > 0;
negative = (Qf < 0)*-1;
Qf = positive + negative;

for i = 1:size(Q)
    if abs(Q(i,1)) < .1
        Q(i,1) = 0;
    end
end


H.Edges.Weight = k;
% HWidths = 5*k/max(k);
% plot(H,'Layout','force','EdgeLabel',Q,'LineWidth',HWidths)


end

