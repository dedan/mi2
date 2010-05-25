
clc
close all
clear all

%% Exercise 5.1 - Toy Data

% load data
load toypca/pca_data.dat
res = struct;


% remove column means (centering)
dat = bsxfun(@minus,pca_data,mean(pca_data,1));

% plot the datapoints
res(1).dat      = dat;
res(1).title    = 'original data';

% compute principal components
[pcs, lambdas]  = eig(cov(dat));
[throw, idx]    = sort(diag(lambdas),'descend');
pcs             = pcs(:,idx);

% transform data to pca space
proj            = dat * pcs;
res(2).dat      = proj;
res(2).title    = 'data in pc space';

% backprojection (reconstruction) from pc space
res(3).dat      = proj * pcs';
res(3).title    = 'backprojected data (full)';

% reconstrucdtion from single pcs
res(4).dat      = proj(:,1) * pcs(:,1)';
res(4).title    = 'backprojected data (pc 1)';

res(5).dat      = proj(:,2) * pcs(:,2)';
res(5).title    = 'backprojected data (pc 2)';

% ploting
for i = 1:length(res)
    subplot(length(res),1,i)
    scatter(res(i).dat(:,1), res(i).dat(:,2));
    title(res(i).title);
    if i~=2
        axis([min(dat(:,1)) max(dat(:,1)) min(dat(:,2)) max(dat(:,2))])
    end
end


%% Exercise 5.2 Image Data

patch_size = 8;                                                 % patch size
num_patches = 500;                                              % number of patches to generate

nat_patches = get_random_patches('n',num_patches,8)';           % extract patches 
b_patches   = get_random_patches('b',num_patches,8)';            

% this time use princomp as centering is already build in
nat_pcs     = princomp(nat_patches);
b_pcs       = princomp(b_patches);

figure(2)
title('pcs of natural images');
for i =1:20
    subplot(5,5,i);
    imagesc(reshape(nat_pcs(:,i),patch_size, patch_size));
    colormap gray
    axis off
end

figure(3)
title('pcs of images of buildings');
for i =1:20
    subplot(5,5,i);
    imagesc(reshape(b_pcs(:,i),patch_size, patch_size));
    colormap gray
    axis off
end

% the pcs of the images of buildings show more rectangular angle
% orientation

