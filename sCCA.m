function Wa = sCCA(X,Y,nIter,comp,thr) 
    q = size(X,2);    
    p = size(Y,2);
    X = normalize(X',comp);    X = X';
    Y = normalize(Y',comp);    Y = Y';
    c1 = thr/sqrt(q);
    c2 = thr/sqrt(p);

%     XtX_reg = X' * X + 1e-8 * eye(q);   % Regularized X'X
%     YtY_reg = Y' * Y + 1e-8 * eye(p); % Regularized Y'Y
%     XtY = X' * Y;                        % Cross-covariance X'Y
% 
%     % Eigenvalue decomposition for X'X + lambda I
%     [U_X, Lambda_X] = eig(XtX_reg, 'vector'); % Vector form for efficiency
%     Lambda_X_inv_sqrt = diag(1 ./ sqrt(Lambda_X)); % Inverse square root of eigenvalues
%     XtX_inv_sqrt = U_X * Lambda_X_inv_sqrt * U_X'; % Inverse square root matrix
% 
%     % Eigenvalue decomposition for Y'Y + lambda I
%     [U_Y, Lambda_Y] = eig(YtY_reg, 'vector'); % Vector form for efficiency
%     Lambda_Y_inv_sqrt = diag(1 ./ sqrt(Lambda_Y)); % Inverse square root of eigenvalues
%     YtY_inv_sqrt = U_Y * Lambda_Y_inv_sqrt * U_Y'; % Inverse square root matrix
% 
%     K = XtX_inv_sqrt * XtY * YtY_inv_sqrt;

    K=X'*Y;
    for r = 1:comp
        u = ((1/p)*ones(p,1))/norm((1/p)*ones(p,1));
        v = ((1/q)*ones(q,1))/norm((1/q)*ones(q,1));
        for k=1:nIter
            vl = v;  
            v = sign(K*u).*max(0, bsxfun(@minus,abs(K*u),c1/2));
            u = sign(K'*v).*max(0, bsxfun(@minus,abs(K'*v),c2/2));
            if sum(v)~=0
                v = v/norm(v);
                u = u/norm(u);
            end
            if abs(norm(v - vl,2))<1e-3 
                break;
            end
        end
        K = K - (v'*K*u)*v*u';
        Wa(:,r) = v;
    end
end



function X = normalize(X,K)
% Each row of X is an attribute. 
% The columns of X are samples
X = X - repmat(mean(X,2),1,size(X,2));
[Ux,~,Vx] = svds(X,K);
X = Ux*Vx';
end
