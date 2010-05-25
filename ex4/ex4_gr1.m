

load toypca/pca_data.dat

% plot the datapoints
subplot 221
scatter(pca_data(:,1), pca_data(:,2)); 
title('original data')

% compute principal components
dat = bsxfun(@minus,pca_data,mean(pca_data,1));
[pcs, lambdas]  = eig(cov(dat));
[throw, idx]    = sort(diag(lambdas),'descend');
pcs             = pcs(:,idx);

% transform data to pca space
proj = dat * pcs;

% datapoints in pc coordinate system
subplot 222
scatter(proj(:,1), proj(:,2))

subplot 223
[t bla] = princomp(pca_data);
scatter(bla(:,1), bla(:,2))

