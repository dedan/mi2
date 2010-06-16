
close all
clear

%% part 1 init

% load the sounds
load sounds/sound1.dat
load sounds/sound2.dat
s = [sound1'; sound2'];

% settings 
N   = 2;
eta = 0.1;
magic_number = 0.001;

% create mixing matrix and mix sources
A = rand(N);
x = A*s;
p = size(x,2);


% permute columns of data matrix
x_perm = x(:,randperm(p));


% correlation coefficient
c = corrcoef(s, x);
disp(['correlation signal vs. sources: ' num2str(c(2,1))]);
c = corrcoef(s, x_perm);
disp(['correlation signal (permuted) vs. sources: ' num2str(c(2,1))]);


% center the data
x = bsxfun(@minus, x, mean(x,2));

% initialize unmixing matrix w randomly
w = rand(N);


%% part 2 optimization

% transform data by logistic function
f = @(x) 1./(1+exp(-x));

changed = true;
while changed
   
    sam = zeros(N);
    for i = 1:p
        y   = w*x_perm(:,i);
        sam = sam + (1-2*f(y)) * y';
    end
    sam = sam ./ p;
    d_w = eta * ((eye(N)+sam) * w);
    w   = w + d_w;
    if norm(d_w) < magic_number
        changed = false;
    end
end


%% results
% Play and plot the original sounds, the mixed sources (before and after permutation) and the recovered signals.
% Calculate the correlations (as above) between the true sources and the estimations.

% recover the signal with the unmixing matrix learned on the permuted mixed
% signals
s_recov = w*x;

% correlation coefficient
c = corrcoef(s, s_recov);
disp(['correlation sources vs. rec. sources: ' num2str(c(2,1))]);

%figure(1)
subplot 611
plot(s(1,:))
soundsc(s(1,:));
subplot 612
plot(s(2,:))
soundsc(s(2,:));

