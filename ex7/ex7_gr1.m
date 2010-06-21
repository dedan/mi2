
clear
close all

load distrib.mat

n_bins = 100;


A           = [4 3; 2 1];               % mixing matrix

res          = struct;
res(1).s     = laplacian;
res(1).label = 'laplacian';
res(1).color = 'g';
res(2).s     = normal;
res(2).label = 'normal';
res(2).color = 'b';
res(3).s     = uniform;
res(3).label = 'uniform';
res(3).color = 'r';

%% 1


% pca spheering of data 
for i = 1:length(res)
    res(i).x         = A * res(i).s;                                  % mix sources
    res(i).x_center  = bsxfun(@minus, res(i).x, mean(res(i).x,2));     % centering
    
    [~, proj, lam]   = princomp(res(i).x');                           % decorrelation
    res(i).x_decor   = proj';
    
    lam              = 1./sqrt(diag(lam));                           % spheering
    lam(lam==Inf)    = 0;
    res(i).x_spheer  = lam * proj';              
end


% rotation and kurtosis
thetas = 0:pi/50:2*pi;
for i = 1:length(res)
    
    res(i).min_value = Inf;
    res(i).max_value = -Inf;   
    res(i).kurt      = ones(length(thetas),2);
    for j = 1:length(thetas)
        
        rot_matrix          = [cos(thetas(j)) -sin(thetas(j)); ...
            sin(thetas(j)) cos(thetas(j))];
        
        res(i).kurt(j,:)    = kurt((rot_matrix * res(i).x_spheer)');
        
        if res(i).kurt(j,1) < res(i).min_value
            res(i).min_value  = res(i).kurt(j,1);
            res(i).min_rot    = rot_matrix;
        end
        if res(i).kurt(j,1) > res(i).max_value
            res(i).max_value  = res(i).kurt(j,1);
            res(i).max_rot    = rot_matrix;
        end
    end
end


% plots
% Plot the original dataset (sources) and the mixed dataset after the 
% steps 1, 2, 3, 4, and 6 as a scatter plot and display the respective marginal histograms. 
% After step 5 plot the kurtosis value as a function of angle for each dimension.
% Compare the histograms for ?min and ?max for the different distributions.
% 1 mixing 2 centering 3 decorelation 4 spheering 6 min max kurt
for i = 1:length(res)

    figure;
    scatterhist(res(i).s(1,:), res(i).s(2,:));
    title(['original data (' res(i).label ')']);
    
    figure;
    scatterhist(res(i).x(1,:), res(i).x(2,:));
    title(['mixed data (' res(i).label ')']);
    
    figure;
    scatterhist(res(i).x_center(1,:), res(i).x_center(2,:));
    title(['centered data (' res(i).label ')']);

    figure;
    scatterhist(res(i).x_decor(1,:), res(i).x_decor(2,:));
    title(['decorrelated data (' res(i).label ')']);
    
    figure;
    scatterhist(res(i).x_spheer(1,:), res(i).x_spheer(2,:));
    title(['spheered data (' res(i).label ')']);
    
    
    figure;
    min_rot = res(i).min_rot * res(i).x_spheer;
    scatterhist(min_rot(1,:), min_rot(2,:));
    title(['min kurt (' res(i).label ')']);

    figure;
    max_rot = res(i).max_rot * res(i).x_spheer;
    scatterhist(max_rot(1,:), max_rot(2,:));
    title(['max kurt (' res(i).label ')']);
    
    figure;
    hold on
    plot(thetas, res(i).kurt(:,1), 'g');
    plot(thetas, res(i).kurt(:,2), 'b');
    legend({'dim 1', 'dim 2'});
    hold off
    
    
    figure(50);
    
    subplot(4,1,1);
    hold on
    range = min(min_rot(1,:)):(max(min_rot(1,:))-min(min_rot))/n_bins:max(min_rot(1,:));
    bar(range, histc(min_rot(1,:),range), 'FaceColor', res(i).color);
    title('rot min dim 1');
    legend({res.label});
    hold off
    
    subplot(4,1,2);
    hold on
    range = min(min_rot(2,:)):(max(min_rot(2,:))-min(min_rot))/n_bins:max(min_rot(2,:));
    bar(range, histc(min_rot(1,:),range), 'FaceColor', res(i).color);
    title('rot min dim 2');
    legend({res.label});
    hold off

    subplot(4,1,3);
    hold on
    range = min(max_rot(1,:)):(max(max_rot(1,:))-min(max_rot))/n_bins:max(max_rot(1,:));
    bar(range, histc(max_rot(1,:),range), 'FaceColor', res(i).color);
    title('rot max dim 1');
    legend({res.label});
    hold off

    subplot(4,1,4);
    hold on
    range = min(max_rot(2,:)):(max(max_rot(2,:))-min(max_rot))/n_bins:max(max_rot(2,:));
    bar(range, histc(max_rot(1,:),range), 'FaceColor', res(i).color);
    title('rot max dim 2');
    legend({res.label});
    hold off
    
end


%% 2

t       = 0:0.05:50;
s       = NaN(3,length(t));
s(1,:)  = 4*sin(t-3);
s(2,:)  = mod(t+5,10);
s(3,:)  = -14 * cos(2*t);
s(3,1)  = 0;

A = [2 -3 -4; 7 5 1; -4 7 5];

x = A*s;

figure(1);
subplot 211
plot(s');
subplot 212
plot(x');

figure(2);
subplot 211
whitesig = fastica(x, 'only', 'white');
plot(whitesig');

subplot 212
sep = fastica(x);
plot(sep');


%% 3

n_patches   = 2000;
patch_size  = 4;
n           = get_random_patches('n',n_patches,patch_size);

[A W] = fastica(n, 'approach', 'symm', 'g', 'tanh');
% figure;
% for i = 1:20
%     subplot(4,5,i);
%     reshape(
%     
    
    
    