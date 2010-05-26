%clc
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


%% kernel pca


toy_data    = [];
n_data      = 30;                           % datapoints per cluster
sig         = 0.1;                          % sigma of dataclusters
means       = [-0.5 -0.2; 0 0.6; 0.5 0];    % means of clusters
n_eig       = 16;                           % number of eigenvectors

% create toy data
for i = 1:size(means,1)
   toy_data = [toy_data randn(2,n_data)*sig + repmat(means(i,:),n_data,1)']; %#ok<AGROW>
end

% compute kernel matrix
p = length(toy_data);
kernel_matrix = zeros(p);
for i = 1:p
    for j = 1:p
        kernel_matrix(i,j) = exp(-norm(toy_data(:,i)-toy_data(:,j))^2/(2*sig^2));
    end
end

% normalize kernel matrix
j = ones(p);
k = kernel_matrix - 1/p * j * kernel_matrix - 1/p * kernel_matrix * j + 1/p^2 * j * kernel_matrix * j;

% compute the eigenvectors of k and sort them
[a_tilde,lambdas]   = eig(k);
[lambdas, idx]      = sort(diag(lambdas),'descend');
a_tilde             = a_tilde(:,idx);

% normalize the eigenvectors
a = NaN(size(a_tilde));
for i = 1:p
    a(:,i) = a_tilde(:,i) / (sqrt(lambdas(i)*p)* norm(a_tilde(:,i)));
end

% creating result struct
range = -1:0.1:1;
res = struct;
for eigk = 1:n_eig
    res(eigk).bla = zeros(length(range));
end

% computing projections for each eigenvector
for i = 1:length(range)
    for j = 1:length(range)
        x = [range(i); range(j)];
        
        % loop over all datapoints
        for beta = 1:p
            % loop over eigenvectors
            for eigk = 1:n_eig
                res(eigk).bla(i,j) = res(eigk).bla(i,j) + norm(a(beta,eigk)) * exp(-norm(toy_data(:,beta)-x)^2/(2*sig^2));
            end
        end
    end
end

% plotting
figure
annotation('textbox','String',{'Kernel-PCS in input space'},'Position',[0.4 0.99 0.5 0],'LineStyle','none');

[x y] = meshgrid(-1:0.1:1, -1:0.1:1);
for i = 1:n_eig
    subplot(4,4,i);
    contour(x, y, res(i).bla)
    hold on
    plot(toy_data(1,:),toy_data(2,:), '.k')
    hold off
end