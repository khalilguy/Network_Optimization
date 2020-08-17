s = [1 1 2 3 3 4 6 7 8 9 4 10 11 12 5];
t = [2 6 3 4 7 5 5 8 9 10 10 11 12 13 13];

nedges= size(s,2);
G = digraph(s,t,1:nedges);
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
figure(1)
plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',LWidths);

nedges = numedges(G)
nnodes = numnodes(G)

it = incidence(G)*-1;
i= (it)';
size(it)
k = G.Edges.Weight;
K = diag(k);

F = sparse(nnodes,1);
F(1,1) = 100;
F(nnodes)= 0;
R = it*K*i;
size(R)
p = 13;
randpressure = sparse(1,nnodes);
randpressure(1,p) = 1;
R(p,:) = randpressure;
P = R\F;




Q = K*i*P;

%Reassigning Conductances
gamma = 1/2;
lambda = 0;
sumK = 0;


%Calculate the conductance constraint
for i = 1:length(k)
    sumK = sumK + k(i,1)^gamma;
end

bigK = sumK^(1/gamma);

%Calculating lambda
for i = 1:length(Q)
    lambda = lambda + (Q(i,1))^((2*gamma)/(gamma + 1));
end

lambda = ((lambda/bigK^gamma)^((gamma + 1)/gamma))/gamma;

newK = zeros(1,length(Q));

 for i = 1:length(Q)
    newK(1,i) = ((Q(i,1)^2)/(lambda*gamma))^(1/(gamma +1));
 end
 
%Recalculating Flows
H = digraph(s,t,newK);
HWidths = 5*newK/max(newK);
figure(2)
plot(H,'Layout','force','EdgeLabel',newK,'LineWidth',HWidths)

it = incidence(H)*-1;
i= (it)';


K = diag(newK);

F = sparse(nnodes,1);
F(1,1) = 100;
F(nnodes)= 0;
R = it*K*i;

p = 13;
randpressure = sparse(1,nnodes);
randpressure(1,p) = 1;
R(p,:) = randpressure;
P = R\F;

newQ = K*i*P;





        

    