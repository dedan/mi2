clear all
clc

im              = imread('testimg.jpg');       % load the image
sigmas          = [0.05 0.1];                  % sigmas for gaussian noise
sample_sizes    = [100 500];                   % sample sizes
min_log_like    = zeros(length(sigmas),length(sample_sizes));       % matrix for min likelihood
h_min_log_like  = zeros(length(sigmas),length(sample_sizes));       % h associated with min likelihood


% iterate over sigmas
for i_sigma = 1:length(sigmas)
    
    % current sigma
    sigma = sigmas(i_sigma);
    
    % iterate over sample sizes
    for i_sample_size = 1:length(sample_sizes)
        
        % current sample size
        sample_size = sample_sizes(i_sample_size);
        
        noise           = randn(size(im)) * sigma * 255;    % generate noise
        noisy_im        = round(double(im)+noise);          % add noise to the image
        noisy_im_vect   = noisy_im(:);

        perms           = randperm(length(noisy_im_vect));      % compute a random permutation of the samples
        sample_indexes  = perms(1:sample_size);                 % sample indexes 
        val_indexes     = perms(sample_size:end);               % validation indexes
        sample_set      = noisy_im_vect(sample_indexes);        % sample set
        val_set         = noisy_im_vect(val_indexes);           % validation set

        figure();

        % show noisy image
        subplot(2,2,1)
        imshow(noisy_im,[0 255]);                  
        title(sprintf('noisy image, sigma = %.3f',sigma));
        
        % show noisy image
        subplot(2,2,2)
        imshow(reshape(sample_set,10,sample_size/10),[0 255]);                  
        title(sprintf('sample set, P = %d',sample_size));
        
        subplot(223)
        hold all

        h_count = 1;                    % counter for looping over h_range
        h_range = [1 10 20 30 40 50];   % h values considered
        labels = {length(h_range)};     % cell array for labels
        log_likelihood = zeros(length(h_range),1);  % log-likelihood vector
        
        % loop over all h values
        for i_h = 1:length(h_range)
            
           h = h_range(i_h);

           x = -h:(255+h);              % x value for the kernel function
           f = zeros(length(x),1);      % probabilty density function

           % loop over points in the patch
           for i=1:length(sample_set)
               xi = sample_set(i);                                  % point in the patch
               k = (1 / sqrt(2*pi) ) * exp(-(x-xi).^2 / (2*h^2));   % kernel evaluated for this point
               f = f+k';                                            % sum up the kernel for this point
           end
           f = f / (length(sample_set) * h);  % final normalization of the pdf

           % compute the log-likelihood of the remaining data
           hist_val_set = histc(val_set,x);                        % count values for each sample of the validation set
           log_likelihood(i_h) = -sum(log(f).*hist_val_set);       % log likelihood

           plot(x,f);                               % plot the pdf
           labels{i_h} =sprintf('h=%d',h);      

        end
        legend(labels)
        title('pdf estimation using gaussian kernel with sigma h')
        xlabel('gray-scale value')
        ylabel('prob.')

        subplot(224)
        plot(h_range,log_likelihood)
        xlabel('h')
        ylabel('- log-likelihood')
        title('- log-likelihood')
        
        % save minimum value of log-likelihood and the corresponding best h
        [m,c] = min(log_likelihood);
        min_log_like(i_sigma,i_sample_size) = m;
        h_min_log_like(i_sigma,i_sample_size) = h_range(c);
        
    end
end

fig = figure();
fig_scale = 0.7;
fig_rel_h = 0.5;
set(fig,'Units','normalized','Position',[(1-fig_scale)/2, (1-fig_scale*fig_rel_h)/2,fig_scale,fig_scale*fig_rel_h])
set(fig,'PaperPositionMode','auto')

% plot minimum log-likelihood
subplot(121)
imagesc(sigmas,sample_sizes,min_log_like);
axis('square')
xlabel('sigma')
ylabel('P')
title('minimum of neg. log-likelihood');
colormap('jet')
set(gca,'xtick',sigmas);
set(gca,'ytick',sample_sizes);
colorbar;

% plot h associated to minimum log-likelihood
subplot(122)
imagesc(sigmas,sample_sizes,h_min_log_like,[1 50]);
axis('square')
xlabel('sigma')
ylabel('P')
title('h with minumum neg. log-likelihood');
colormap('jet')
set(gca,'xtick',sigmas);
set(gca,'ytick',sample_sizes);
colorbar;