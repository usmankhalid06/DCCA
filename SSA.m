function RC = SSA(X,M)
N = size(X,1);
covX = xcorr(X,M-1,'unbiased');
Ctoep=toeplitz(covX(M:end));

Y=zeros(N-M+1,M);
for m=1:M
  Y(:,m) = X((1:N-M+1)+m-1);
end;
Cemb=Y'*Y / (N-M+1);
C = Cemb; 
[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors
PC = Y*RHO;

RC=zeros(N,M);
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;



