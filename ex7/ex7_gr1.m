clc
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

%% Exercise 1 - Kurtosis of Toy Data 


% pca spheering of data 
for i = 1:length(res)
    res(i).x         = A * res(i).s;                                % mix sources
    res(i).x_center  = bsxfun(@minus, res(i).x, mean(res(i).x,2));  % centering
    
    [del, proj, lam]   = princomp(res(i).x');                       % decorrelation
    res(i).x_decor   = proj';
    
    lam              = 1./sqrt(diag(lam));                          % spheering
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
        
        % rotation matrix
        rot_matrix          = [cos(thetas(j)) -sin(thetas(j)); ...
                               sin(thetas(j)) cos(thetas(j))];
                           
        % compute kurtosis of rotated data
        res(i).kurt(j,:)    = kurtosis((rot_matrix * res(i).x_spheer)');
        
        % save minimum kurtosis
        if res(i).kurt(j,1) < res(i).min_value
            res(i).min_value  = res(i).kurt(j,1);
            res(i).min_rot    = rot_matrix;
        end
        
        % save maximum kurtosis
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
% Compare the histograms for kurt_min and kurt_max for the different distributions.
% 1 mixing 2 centering 3 decorelation 4 spheering 6 min max kurt
for i = 1:length(res)

    figure;
    subplot(2,3,1);
    scatterhist(res(i).s(1,:), res(i).s(2,:));
    title(['original data (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_original.eps'));
    
    figure;
    scatterhist(res(i).x(1,:), res(i).x(2,:));
    title(['mixed data (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_mixed.eps'));
    
    figure;
    scatterhist(res(i).x_center(1,:), res(i).x_center(2,:));
    title(['centered data (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_centered.eps'));

    figure;
    scatterhist(res(i).x_decor(1,:), res(i).x_decor(2,:));
    title(['decorrelated data (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_decorellated.eps'));
    
    figure;   
    scatterhist(res(i).x_spheer(1,:), res(i).x_spheer(2,:));
    title(['spheered data (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_spheered.eps'));
    
    figure;
    subplot(1,2,1);
    min_rot = res(i).min_rot * res(i).x_spheer;
    scatterhist(min_rot(1,:), min_rot(2,:));
    title(['min kurt (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_minkurt.eps'));
    
    figure;
    max_rot = res(i).max_rot * res(i).x_spheer;
    scatterhist(max_rot(1,:), max_rot(2,:));
    title(['max kurt (' res(i).label ')']);
    print('-depsc',strcat(res(i).label,'_maxkurt.eps'));
    
    figure;
    hold on
    plot(thetas, res(i).kurt(:,1), 'g');
    plot(thetas, res(i).kurt(:,2), 'b');
    legend('dim 1', 'dim 2');
    print('-depsc',strcat(res(i).label,'_kurtrot.eps'));
    hold off
    
    scrsz = get(0,'ScreenSize');
    figure(50);
    set(gcf,'Position',[1 1 scrsz(3)*0.8 scrsz(4)/2])
    
    subplot(2,2,1);
    hold on
    range = min(min_rot(1,:)):(max(min_rot(1,:))-min(min_rot))/n_bins:max(min_rot(1,:));
    bar(range, histc(min_rot(1,:),range),1, 'FaceColor', res(i).color,'EdgeColor','none');
    xlim([-6 6])
    ylim([0 600])
    title('rot min (dim 1)');
    legend({res.label});
    hold off
    
    subplot(2,2,2);
    hold on
    range = min(min_rot(2,:)):(max(min_rot(2,:))-min(min_rot))/n_bins:max(min_rot(2,:));
    bar(range, histc(min_rot(1,:),range),1, 'FaceColor', res(i).color,'EdgeColor','none');
    xlim([-6 6])
    ylim([0 600])
    title('rot min (dim 2)');
    legend({res.label});
    hold off

    subplot(2,2,3);
    hold on
    range = min(max_rot(1,:)):(max(max_rot(1,:))-min(max_rot))/n_bins:max(max_rot(1,:));
    bar(range, histc(max_rot(1,:),range), 1,'FaceColor', res(i).color,'EdgeColor','none');
    xlim([-6 6])
    ylim([0 600])
    title('rot max (dim 1)');
    legend({res.label});
    hold off

    subplot(2,2,4);
    hold on
    range = min(max_rot(2,:)):(max(max_rot(2,:))-min(max_rot))/n_bins:max(max_rot(2,:));
    bar(range, histc(max_rot(1,:),range), 1,'FaceColor', res(i).color,'EdgeColor','none');
    xlim([-6 6])
    ylim([0 600])
    title('rot max (dim 2)');
    legend({res.label});
    hold off
    print('-depsc','kurthist.eps');
    
end


%% Exercise 2 - Toy Signal Separation 

t = 0:0.05:50; % time

% sources
s       = NaN(3,length(t));
s(1,:)  = 4*sin(t-3);
s(2,:)  = mod(t+5,10);
s(3,:)  = -14 * (cos(2*t)>0);

A = [2 -3 -4; 7 5 1; -4 7 5];   % mixing matrix
x = A*s;                        % mix sources

figure;
subplot 221
plot(t,s');
xlim([0 50]);
title('original sources');

subplot 222
plot(t,x');
xlim([0 50]);
title('mixed signal');

subplot 223
whitesig = fastica(x, 'only', 'white');
plot(t,whitesig');
xlim([0 50]);
title('whitened signal');

subplot 224
sep = fastica(x);
plot(t,sep');
xlim([0 50]);
title('reconstructed sources')

print('-depsc','toysep.eps');

%% Exercise 3 - ICA on Image Patches 

n_patches   = 20000;    % number of patches per category
patch_size  = 16;       % patch size 16x16

% image categories
types = ['n','b','t'];
names = {'nature','building','text'};

% repeat for each category
for i=1:length(types)
    
    n = get_random_patches(types(i),n_patches,patch_size);
    [A W] = fastica(n, 'approach', 'symm', 'g', 'tanh','lastEig',20);
    
    % plot first 20 independent features
    figure;
    annotation('textbox','String',names{i},'Position',[0.45 0.99 0.5 0],'LineStyle','none');
    for j = 1:20
        subplot(4,5,j);
        imagesc(reshape(A(:,j),16,16));
        set(gca,'xTick',[]);
        set(gca,'yTick',[]);
    end
    print('-depsc',strcat(names{i},'_indfeats.eps'));
end
     
    
    
    