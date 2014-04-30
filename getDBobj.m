function [objective, DB] = getDBobj(X,label,m)
%% compute the Davies-Bouldin index
% X is the data represented as rows
% label is the cluster label
% m is the centroids 
k = size(m,1);
S = zeros(1,k);
M = zeros(k);
R = zeros(k,k);

centroid = m;
objective = 0; 
for i = 1:k
    dat = X(label == i,:);
    T = size(dat,1);
    
    tmp = 0;
   for j = 1:T
        tmp = tmp + norm(dat(j,:) - centroid(i,:))^2;
    end
    objective=objective+tmp;
    tmp = (tmp/T)^.5;
    S(i) = tmp;
end

objective = objective/size(X,1);
    
for i = 1:k
    for j = 1:k
        M(i,j) = norm(centroid(i,:) - centroid(j,:));
    end
end


for i = 1:k
    for j = 1:k
        if i == j
            continue;
        end
        R(i,j) = (S(i) + S(j))/M(i,j);
    end
end

D = max(R,[],2);
DB = mean(D);
    
end
