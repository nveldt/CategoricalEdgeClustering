%% First: Install and compile graclus 1.2
%    http://www.cs.utexas.edu/users/dml/Software/graclus.html
addpath('graclus1.2/matlab/')

%% Run Graclus
load A_sym
n = size(A,1);
ks = 2:40;
Cmets = zeros(n,numel(ks));

for t = 1:numel(ks)
    k = ks(t);
    [c, obj] = graclus(A, k);
    Cmets(:,t) = c;
end

save('Graclus_clusterings_2to40.mat','Cmets')
