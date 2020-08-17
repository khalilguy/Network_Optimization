function [C, s, t, location, circle] = create_circle(new_nodes, num_nodes, radius, tol)
%put a source at the middle of the graph
gridmiddle = [sum(new_nodes(:,1))/num_nodes, sum(new_nodes(:,2))/num_nodes];
for i = 1:num_nodes
    diffs(i, 1) = sum((new_nodes(i, :) - gridmiddle).^2);
end
[val, idx] = min(diffs);
source = new_nodes(idx, :);

dist_center = zeros(num_nodes, 1);
for i = 1:num_nodes
    dist_center(i, 1) = sqrt((new_nodes(i, 1)-source(1,1)).^2+(new_nodes(i, 2)-source(1,2)).^2);
end

in_circle = dist_center< radius + tol;
indices_in = find(in_circle);

circle = new_nodes(indices_in,:);

[c_nodes, ncol] = size(circle);
[x_c,y_c] = meshgrid(circle(:,1),circle(:,2));
dist_c = sqrt((x_c-x_c.').^2+(y_c-y_c.').^2);

dist_c(triu(true(c_nodes))) = NaN;
Cis_adjac = dist_c< 1+ tol;
[t, s] = find(Cis_adjac);

%%
% hold on
% scatter(circle(:,1), circle(:,2))
% th = 0:pi/50:2*pi;
% xunit = radius * cos(th) + source(1,1);
% yunit = radius * sin(th) + source(1,2);
% h = plot(xunit, yunit);
% hold off

%%
C = digraph(s, t);
% plot(digraph(s,t), 'Layout', 'force')
circle_source = indices_in == idx;
location = find(circle_source);