clear;
close all; 
clc;



%% parameters
sp = simtb_create_sP('exp_params_aod');
N = sp.nT; %number of time samples 
nV = sp.nV; %sqrt of number of voxels
nSRCS = sp.nC; %number of sources
nIter = 30; %algorithm iterations
tstd  = sqrt(0.6); %0.6 is the varaince
sstd  = sqrt(0.015);
Dp = dctbases(N,N); %dct basis dictionary
nf = sqrt(70*15.49); 

% for j=1:1
j = 1;    

    sv = 10;
    sp.SM_spread =sv+0.05*randn(N,nV*nV); %0.0001
    SM = simtb_makeSM(sp,1);   % Create spatial maps
    TC = zscore(simtb_makeTC(sp,1));  % Create TCs 
    rng('default'); 
    rng('shuffle') % random number generator
    Y= (TC+tstd*randn(N,nSRCS))*(SM+sstd*randn(nSRCS,nV^2));
    Y= Y-repmat(mean(Y),size(Y,1),1);    



    %% PCA+PMD
    tic
    K1 = 8; K2 = 16;
    [Y2,~,~]=svds(Y,K1);
    Y2 = Y2 *diag(1./sqrt(sum(Y2 .*Y2)));
    lambda1 = 0.5;  
    [Wx1,~,~] = pmd_rankK(Y',Y2',K2,lambda1);
    X{1} = Wx1';
    D{1} = Y*X{1}';
    toc

    %% SPCA+PMD
    tic;
    K1 = 8; K2 = 16;
%     [tmpX,~] = spca(Y,[],K1,Inf,-150, 1000, 1e-6); % 2000
    lambda =0.5; gamma=lambda*ones(1,K1);
    tmpX=GPower(Y,gamma,K1,'l1',0);
    Y2 = Y*tmpX;
    Y2 = Y2*diag(1./sqrt(sum(Y2.*Y2)));
    lambda2 = 0.5;
    [Wx1,~,~] = pmd_rankK(Y',Y2',K2,lambda2);
    X{2} = Wx1';
    D{2} = Y*X{2}';
    toc

    %% ICA+PMD
    tic;
    K1 = 8;  K2 = 16;
    [G,~,~] = svds(Y,K1);
    Ss = G'*Y;
    [SSs,A,~] = fastica(Ss, 'numOfIC', K1,'approach','symm', 'g', 'tanh','verbose', 'off');
    Y2 = Y*SSs'; 
    Y2 = Y2 *diag(1./sqrt(sum(Y2 .*Y2)));
    lambda2 = 0.5; 
    [Wx1,~,~] = pmd_rankK(Y',Y2',K2,lambda2);
    X{3} = Wx1';
    D{3} = Y*X{3}';
    toc

    %% sICA+PMD
    tic;
    K1 = 8;  K2 = 16;
    [G,~,~] = svds(Y,K1);
    Ss = G'*Y;
    [WW,~,~,~,~] = ICA_EBM_Sparse(Ss,5,10^6);
    Y2 = Y*(WW*Ss)'; 
    Y2 = Y2 *diag(1./sqrt(sum(Y2 .*Y2)));
    lambda2 = 0.5; 
    [Wx1,~,~] = pmd_rankK(Y',Y2',K2,lambda2);
    X{4} = Wx1';
    D{4} = Y*X{4}';
    toc;

    %% DCCAPMD
    tic
    K = 16;
    [Y2,~]= kmeans_clustering(Y,8,12,32,1); %16
    Y2 = Y2 *diag(1./sqrt(sum(Y2 .*Y2)));
    lambda2 = 0.5; 
    [Wx2,~]= pmd_rankK(Y',Y2',K,lambda2);
    X{5} = Wx2';
    D{5} = Y*X{5}';    
    for i =1:size(X{5},1); X{5}(i,:) = spatial(abs(X{5}(i,:)),nV,nV); end 
    for i =1:size(D{5},2)
        U4 = SSA(D{5}(:,i),50);
        AA = zeros(size(U4,2),1);
        [~,bb]= sort(abs(U4'*D{5}(:,i)),'descend');
        ind = bb(1:18);
        AA(ind,1)= (U4(:,ind)'*U4(:,ind))\U4(:,ind)'*D{5}(:,i);
        AA = AA./norm(U4*AA);
        D{5}(:,i) = U4*AA;
    end    
    toc

    warning('off','all')
    %% DCCAPCA
    tic;
    K = 16; %nnz(X{6})/(size(X{6},1)*size(X{6},2))
    spa1 = (36/nf)*nf; %40percent at 8
    spa2 = 0.5;%2-
    spa3 = 4;%8
    [Dc,Xc]= kmeans_clustering(Y,8,12,32,1); %18
    Dc2 = filter([0 1],1,Dc);
    for i =1:size(Xc,1); Xc2(i,:) = spatial(abs(Xc(i,:)),nV,nV); end 
    Xc3 = [Xc; Xc2]; %Xc2
    Dc3 = [Dc Dc2];  %Dc2
    iDc=[];
    for i =1:size(Dc3,2)
        tmpiDc = SSA(Dc3(:,i),50);
        iDc = [iDc tmpiDc(:,1:5)]; %5
    end
    iDc = iDc*diag(1./sqrt(sum(iDc.*iDc)));
    [D{6},X{6},Err,A2,B2,EE]= my_DPCA_sim_3(Y,iDc,Xc3,K,spa1,spa2,spa3,nIter);
    for i =1:size(D{6},2)
        U4 = SSA(D{6}(:,i),50);
        AA = zeros(size(U4,2),1);
        [~,bb]= sort(abs(U4'*D{6}(:,i)),'descend');
        ind = bb(1:20);
        AA(ind,1)= (U4(:,ind)'*U4(:,ind))\U4(:,ind)'*D{6}(:,i);
        AA = AA./norm(U4*AA);
        D{6}(:,i) = U4*AA;
    end  

    %%
    sD{1} = TC; 
    sX{1} = SM;
    nA = 7;
    for jj =1:nA-1
%         [sD{jj+1},sX{jj+1},indd{jj}]=sort_TSandSM_temporal(TC,D{jj},X{jj});

        [sD{jj+1},sX{jj+1},indd{jj}]=sort_TSandSM_spatial(TC,SM,D{jj},X{jj},nSRCS);
        for ii =1:nSRCS
            TCcorr(jj+1,ii,j) =abs(corr(TC(:,ii),D{jj}(:,indd{jj}(ii))));
            SMcorr(jj+1,ii,j) =abs(corr(SM(ii,:)',X{jj}(indd{jj}(ii),:)'));
        end
    end

nA = 7;
vec = 1:1;
ccTC = mean(TCcorr(:,:,vec),3);
ccSM = mean(SMcorr(:,:,vec),3);

TT(1,:) = mean(ccTC');
TT(2,:) = mean(ccSM');


TP = zeros(nA,nSRCS); FP = zeros(nA,nSRCS); FN = zeros(nA,nSRCS); F_score = 0; thr = 0.01;
for i =1:nA-1
    for jjj=1:nSRCS
        SM2(jjj,:) =SM(jjj,:)/norm(SM(jjj,:)); 
        [~, inddd(jjj)]  = max(abs(corr(abs(SM2(jjj,:)'),abs(X{i}'))));
        XX{i}(inddd(jjj),:) = X{i}(inddd(jjj),:)/norm(X{i}(inddd(jjj),:));
        TP(i,jjj) = TP(i,jjj) +sum(abs(SM2(jjj,:))>=thr & abs(XX{i}(inddd(jjj),:))>=thr);
        FP(i,jjj) = FP(i,jjj) +sum(abs(SM2(jjj,:))<=thr & abs(XX{i}(inddd(jjj),:))>=thr);
        FN(i,jjj) = FN(i,jjj) +sum(abs(SM2(jjj,:))>=thr & abs(XX{i}(inddd(jjj),:))<=thr); 
        F_score(i,jjj) = (2*(TP(i,jjj)))/(2*(TP(i,jjj))+sum(FP(i,jjj))+(FN(i,jjj)));
    end
end
F_score = [zeros(1,nSRCS); F_score];
mean(F_score)



TT(3,:) =mean(F_score'); 


for k =1:6
    for i =1:8
        aa = reshape(abs(sX{1}(i,:)),nV,nV); aa = aa/norm(aa); aa(aa<=0.00001)=0; aa(aa>0.00001)=1;
        bb = reshape(abs(sX{k+1}(i,:)),nV,nV); bb = bb/norm(bb); bb(bb<=0.00001)=0; bb(bb>0.00001)=1;
        DS(i,k) = dice(aa,bb);
    end
    
end
DS = DS';
DS = [zeros(1,8); DS];
TT(4,:) =mean(DS');
TT    

f = figure; f.Position = [10 250 1650 650];  my_subplots(nA, nSRCS,nV,nV,TCcorr(:,:,j),F_score,sD,sX);





