function [out_node, in_node, sinks,new_nodes,num_nodes] = hexagonal_graph_fixed_sinks(N,M,source,sinks)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tol = 0.01;
x_hexagon=[-1 -0.5 0.5 1 0.5 -0.5 -1];
y_hexagon=[0 -sqrt(3)/2 -sqrt(3)/2 0 sqrt(3)/2 sqrt(3)/2 0];


i = 1;

for nn=0:N
    for mm=0:M
        nodes(i:i+6, 1) = x_hexagon'+3*nn;
        nodes(i:i+6, 2) = y_hexagon'+sqrt(3)*mm;
        i = i + 7;
    end
end


%%

sz = size(nodes);
num_nodes = sz(1, 1);


[x_g,y_g] = meshgrid(nodes(:,1),nodes(:,2));
dist_g = sqrt((x_g-x_g.').^2+(y_g-y_g.').^2);

dist_g(triu(true(num_nodes))) = NaN;
is_repeat = dist_g<tol;

n = 1;

[repeat,og] = find(is_repeat);

szR = size(repeat);
num_repeats = szR(1,1);
for i = 1:num_repeats
    m = repeat(i, 1);
    nodes(m, :) = [100000,100000];
end

for i = 1:num_nodes
    if nodes(i, 1) ~= 100000
        new_nodes(n, :) = nodes(i, :);
        n = n + 1;
    end
end

[num_nodes, ncol] = size(new_nodes);
[x_g,y_g] = meshgrid(new_nodes(:,1),new_nodes(:,2));
dist_g = sqrt((x_g-x_g.').^2+(y_g-y_g.').^2);

dist_g(triu(true(num_nodes))) = NaN;
is_adjac = dist_g< 1+ tol;
[in_node, out_node] = find(is_adjac);




H = digraph(out_node, in_node);

I = incidence(H)';

%Defining sinks and sorces so that the edges will orient themselves in the
%direction of the sink

num_sinks = size(sinks,1);

%Obtaining the x and y coordinates for each of the sinks
sink_positions = zeros(num_sinks,2);
for i = 1 : num_sinks
    sink_positions(i,:) = new_nodes(sinks(i,1),:);
end


num_edges = size(I,1);
num_nodes = size(I,2);

%Calculating the distance that each sink is away from every other node
sink_distances = zeros(num_nodes,num_sinks);

for j = 1 : num_sinks
    xs = sink_positions(j,1);
    ys = sink_positions(j,2);
    for i = 1:num_nodes
        sink_distances(i,j) = sqrt((xs- new_nodes(i,1))^2 + (ys - new_nodes(i,2))^2);
    end
end

%The next step is to determine which sink a node is closest to
% closest_sink_distance = zeros(num_nodes,num_sinks);
% 
% 
% for i = 1:num_nodes
%     distances = sink_distances(i,:);
%     position = 0;
%     for j = 1 : num_sinks
%         smallest = distances(1,1);
%         if distances(1,j) < smallest
%             smallest = distances(1,j);
%             position = j;
%         end
%     end
%     closest_sink_distance(i,position) = smallest;
% end
            
%Swap the which node is a target and a souce node if the source node is closer to the sink than the target node        


for j = 1:num_sinks
    distances_for_one_sink = sink_distances(:,j);
    for i = 1:num_edges
      if distances_for_one_sink(out_node(i,1),1) < distances_for_one_sink(in_node(i,1),1)
         if (distances_for_one_sink(out_node(i,1),1)) ~=0 || (distances_for_one_sink(in_node(i,1),1) ~= 0)
            new_out_node = in_node(i,1);
            in_node(i,1) = out_node(i,1);
            out_node(i,1) = new_out_node;
          end 
      end
    end
end

%We also need to fix the orientation of the edges around the source and
%sinks

for i =1:num_edges
    if in_node(i,1) == source
         new_out_node = in_node(i,1);
         in_node(i,1) = out_node(i,1);
         out_node(i,1) = new_out_node;
    end
    for j = 1 : num_sinks
        if out_node(i,1) == sinks(j,1)
            new_in_node = out_node(i,1);
            out_node(i,1) = in_node(i,1);
            in_node(i,1) = new_in_node;
        end 
    end       
end
end

