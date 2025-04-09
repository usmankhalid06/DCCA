function [D,X,Err,A,B,E]= my_DPCA_sim_3(Y,Dp,Xp,K,spa1,spa2,spa3,nIter)
    K1 = size(Dp,2);
    K2 = size(Xp,1);
    A = eye(K1,K);
    B = zeros(K,K2);
    D = Dp*A;
    X = B*Xp;
    fprintf('Iteration:     ');
    R = Y;
    gridSize = 70;
    minClusterSize = 1000; %- Minimum size of clusters to retain

%         sparseAdj = sparse_knn_adjacency(coords, k);

    for iter=1:nIter
        fprintf('\b\b\b\b\b%5i',iter);
        Do = D;
        for j =1:K
            X(j,:) = 0; A(:,j) = 0; B(j,:) = 0;
            E = R-D*X;
            
            xk = D(:,j)'*E;
            thr1 = spa1./abs(xk); 
            xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr1/2));
            B(j,:) = xkk*Xp'/(Xp*Xp');
            X(j,:) = B(j,:)*Xp;

            X(j,:) = firm_thresholding_nonadaptive(X(j,:), spa2/2, spa3/2);  
            
            rInd = find(X(j,:));
            if (length(rInd)>1)
                tmp3 = E(:,rInd)*X(j,rInd)'; 
                A(:,j)= (Dp'*Dp)\Dp'*tmp3;    
                A(:,j) = A(:,j)./norm(Dp*A(:,j));
                D(:,j) = Dp*A(:,j);
            end          
        end     


%         E = Y - D*X;
%         E= sign(E).*max(0, bsxfun(@minus,abs(E),3/2));
%         R = Y-E;

        Err(iter) = sqrt(trace((D-Do)'*(D-Do)))/sqrt(trace(Do'*Do));

    end
end




% Helper: Remove Small Clusters Using Adjacency
function a = remove_small_clusters_adj(a, sparseAdj, minClusterSize)
    % Remove small clusters from spatial map a based on adjacency
    clusters = conncomp(graph(sparseAdj)); % Identify connected components
    for c = unique(clusters)'
        clusterIdx = find(clusters == c);
        if numel(clusterIdx) < minClusterSize
            a(clusterIdx) = 0; % Remove small clusters
        end
    end
end

