
close all
clear

%% part 1 init

% collect the results in here
res = struct;

% load the sounds
load sounds/sound1.dat
load sounds/sound2.dat
s               = [sound1'; sound2'];
res(1).s        = s;
res(1).label    = 'orignal signal';

% settings 
N               = 2;
eta             = 0.001;                    % learning rate
magic_number    = 0.001;                    % convergence criterion

% create mixing matrix and mix sources
A               = rand(N);
x               = A*s;
res(2).s        = x;
res(2).label    = 'mixed signal';

p               = size(x,2);


% permute columns of data matrix
x_perm          = x(:,randperm(p));
res(3).s        = x_perm;
res(3).label    = 'mixed signal (permuted)';


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

% logistic function
f = @(x) 1./(1+exp(-x));

changed = true;             % stopping criterion
while changed
   
    % iteration over patterns (online)
    for i = 1:p
        y   = w*x_perm(:,i);
        d_w = eta * ((eye(N)+(1-2*f(y)) * y') * w);
        w   = w + d_w;
        if norm(d_w) < magic_number
            changed = false;
        end
    end
end


%% results

% recover the signal with the unmixing matrix learned on the permuted mixed
% signals
s_recov         = w*x;
res(4).s        = s_recov;
res(4).label    = 'recovered signal';


% correlation coefficient
c = corrcoef(s, s_recov);
disp(['correlation sources vs. rec. sources: ' num2str(c(2,1))]);

% the plots
figure(1)
for i = 1:length(res)
    subplot(length(res),2,(i*2)-1)
    plot(res(i).s(1,:))
    soundsc(res(i).s(1,:));
    title(res(i).label);
    
    subplot(length(res),2,i*2)
    plot(res(i).s(2,:))
    soundsc(res(i).s(2,:));
    title(res(i).label);
end
