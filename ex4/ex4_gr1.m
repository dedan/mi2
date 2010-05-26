clc
close all
clear all

%% Exercise 5.1 - Toy Data

% load data
load toypca/pca_data.dat

[v,l] = princomp(pca_data);

% center data
mu = mean(pca_data);
cent_data(:,1) = pca_data(:,1) - mu(1);   
cent_data(:,2) = pca_data(:,2) - mu(2);   

% compute principal components and sort them by relevance
[pcs, lambdas]  = eig(cov(cent_data)); 
[throw, idx]    = sort(diag(lambdas),'descend');
pcs             = pcs(:,idx);

% transform data to pca space
coeff = cent_data * pcs;

fig = figure();
fig_scale = .5;
fig_rel_h = 1; %relative height
set(fig,'Units','normalized','Position',[(1-fig_scale)/2, (1-fig_scale*fig_rel_h)/2,fig_scale,fig_scale*fig_rel_h])
set(fig,'PaperPositionMode','auto')

% plot the datapoints
subplot 121
plot(pca_data(:,1), pca_data(:,2),'.k'); 
axis([0 4 -4 0])
axis square
title('original data')

% datapoints in pc coordinate system
subplot 122
plot(coeff(:,1), coeff(:,2),'.k')
axis square
title('PC coordinate system')

% reconstruct data
rec_data = size(pca_data);      % totally reconstructed data
rec_data_v1 = size(pca_data);   % data reconstructed only with the first pca
rec_data_v2 = size(pca_data);   % data reconstructed only with the second pca

for i=1:size(pca_data,1)
    rec_data(i,:) = coeff(i,:)*pcs + mu;
    rec_data_v1(i,:) = coeff(i,1)*pcs(:,1) + mu';
    rec_data_v2(i,:) = coeff(i,2)*pcs(:,2) + mu';  
end

fig = figure();
fig_scale = .9;
fig_rel_h = .5; %relative height
set(fig,'Units','normalized','Position',[(1-fig_scale)/2, (1-fig_scale*fig_rel_h)/2,fig_scale,fig_scale*fig_rel_h])
set(fig,'PaperPositionMode','auto')

subplot(1,3,1)
plot(rec_data(:,1), rec_data(:,2),'.k'); 
axis([0 4 -4 0])
axis square
title('reconstructed data')

subplot(1,3,2)
plot(rec_data_v1(:,1), rec_data_v1(:,2),'.k'); 
axis([0 4 -4 0])
axis square
title('reconstructed data using v1')

subplot(1,3,3)
plot(rec_data_v2(:,1), rec_data_v2(:,2),'.k'); 
axis([0 4 -4 0])
axis square
title('reconstructed data using v2')


%% Exercise 5.2 Image Data

patch_size = 8;                                                 % patch size
num_patches = 5000;                                             % number of patches to generate

nat_patches = get_random_patches('n',num_patches,8)';           % extract patches 
b_patches   = get_random_patches('b',num_patches,8)';            

% this time use princomp as centering is already build in
nat_pcs     = princomp(nat_patches);
b_pcs       = princomp(b_patches);

min_n = min(min(nat_pcs(:,1:20)));
max_n = max(max(nat_pcs(:,1:20)));

min_b = min(min(b_pcs(:,1:20)));
max_b = max(max(b_pcs(:,1:20)));

figure()
for i =1:20
    subplot(4,5,i);
    imagesc(reshape(nat_pcs(:,i),patch_size, patch_size),[min_n max_n]);
    colormap gray
    axis off
end
annotation('textbox','String',{'PCS of natural images'},'Position',[0.4 0.99 0.5 0],'LineStyle','none');


figure()
for i =1:20
    subplot(4,5,i);
    imagesc(reshape(b_pcs(:,i),patch_size, patch_size),[min_b max_b]);
    colormap gray
    axis off
end
annotation('textbox','String',{'PCS of buildings'},'Position',[0.4 0.99 0.5 0],'LineStyle','none');

% the pcs of the images of buildings show more rectangular angle
% orientation

