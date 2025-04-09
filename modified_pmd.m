function [D,Wx] = modified_pmd(Dq,Dq2,lambda,K,M)


[Wx,~,di] = pmd_rankK(Dq,Dq2,K,lambda);
Wx=Wx';
D = Dq'*Wx';
%D = Dq *diag(1./sqrt(sum(Dq .*Dq)));
for i =1:size(D,2)
    U4 = SSA(D (:,i),M);
    AA = zeros(size(U4,2),1);
    [~,bb]= sort(abs(U4'*D (:,i)),'descend');
    ind = bb(1:18);
    AA(ind,1)= (U4(:,ind)'*U4(:,ind))\U4(:,ind)'*D(:,i);
    AA = AA./norm(U4*AA);
    D(:,i) = U4*AA;
end
