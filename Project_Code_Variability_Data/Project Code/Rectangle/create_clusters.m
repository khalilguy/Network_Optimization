
load particular_configuration2by2.mat distance_matrix ordered_configs ordered_freq_configs
load particular_flow_configuration2by2.mat total_flow_configs_idx total_flow_configs ordered_flow_configs



flow1 = abs(ordered_flow_configs);
flow = sort(flow1,2, 'descend');
hdist = pdist(flow,'euclidean');
hsquare = squareform(hdist);
dist_links = linkage(hsquare,'average');

c = cophenet(dist_links,hdist)
[dendo,T] = dendrogram(dist_links,13);
figure('Name', 'Cluster vs Number of Edges Used');
[sortT,sortTinds] = sort(T,'ascend');
plot(sortT,sum(flow(sortTinds,:) > 1,2));
xlabel('Cluster')
ylabel('Number of Edges Used')

% 'euclidean'	
% Euclidean distance (default).
% 
% 'squaredeuclidean'	
% Squared Euclidean distance. (This option is provided for efficiency only. It does not satisfy the triangle inequality.)
% 
% 'seuclidean'	
% Standardized Euclidean distance. Each coordinate difference between observations is scaled by dividing by the corresponding element of the standard deviation, S = nanstd(X). Use DistParameter to specify another value for S.
% 
% 'mahalanobis'	
% Mahalanobis distance using the sample covariance of X, C = nancov(X). Use DistParameter to specify another value for C, where the matrix C is symmetric and positive definite.
% 
% 'cityblock'	
% City block distance.
% 
% 'minkowski'	
% Minkowski distance. The default exponent is 2. Use DistParameter to specify a different exponent P, where P is a positive scalar value of the exponent.
% 
% 'chebychev'	
% Chebychev distance (maximum coordinate difference).
% 
% 'cosine'	
% One minus the cosine of the included angle between points (treated as vectors).
% 
% 'correlation'	
% One minus the sample correlation between points (treated as sequences of values).
% 
% 'hamming'	
% Hamming distance, which is the percentage of coordinates that differ.
% 
% 'jaccard'	
% One minus the Jaccard coefficient, which is the percentage of nonzero coordinates that differ.
% 
% 'spearman'	
% One minus the sample Spearman's rank correlation between observations (treated as sequences of values).