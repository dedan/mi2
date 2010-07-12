%% Exercise 1

%close all
clear
clc

N = 6;                          % number of nodes
beta = 0.1;                     % inverse of temperature
tau = 1.1;                      % annealing factor
tmax = 100;                      % maximum time


% randomly initialize initial states
S_init = ones(1,N);
S_init(rand(1,N) > 0.5) = -1;

S = S_init;

% randomly initialize weights
W = randn(N,N);
W = W+W';
W = W-diag(diag(W));


betaVect    = zeros(1,tmax);    % to save betas over time
EVect       = zeros(1,tmax);    % to save total network energy over time

% iterate over time
for t=1:tmax
    random_node = randi(N);     % select random node
    
    % compute local energy of this node
    E           = -0.5 *S(random_node) * (W(random_node,:)*S');
    
    % compute probability to change (random_node)
    p           = 1 / (1-exp(-beta*2*E));
    if(rand < p)
        S(random_node) = -S(random_node);
    end
    beta        = tau*beta;     % update beta
    betaVect(t) = beta;         % store results over time
    EVect(t)    = -0.5 * sum(sum(W .* (S' * S)));
end
final_S = S;

figure(1)
subplot 211
plot(EVect)
title('Energy');
xlabel('time');
ylabel('energy');
subplot 212
plot(1./betaVect)
title('beta');
xlabel('time');
ylabel('temperature');


% compute energy for all possible states
E_all = NaN(1,N);
betas = 0.1:0.1:0.5;
P_all = NaN(length(betas), N);
for i=0:2^N-1
    S               = bitget(i, N:-1:1);
    S(S==0)         = -1;
    E_all(i+1)      = -0.5 * sum(sum(W .* (S' * S)));
    P_all(:,i+1)    = exp(-betas * E_all(i+1))';
end

P_all = P_all ./ repmat(sum(P_all,2),1,2^N);

figure(2)
subplot 611;
bar(E_all);
title('Energy of all States');
xlabel('states');
ylabel('energy');

for i=1:length(betas)
    subplot(6,1,i+1);
    bar(P_all(i,:));
    xlabel('states');
    ylabel('probability');
end
    

%% Exercise 2

N = 6;                          % number of nodes
beta = 0.1;                     % inverse of temperature
tau  = 1.1;                    % annealing factor
tmax = 100;                     % maximum time

betaVect    = zeros(1,tmax);    % to save betas over time
EVect       = zeros(1,tmax);    % to save total network energy over time
S = S_init;


% iterate over time
for t=1:tmax
    random_node = randi(N);     % select random node

    e       = W * S';           % compute the mean field
    S       = tanh(beta * e');  % update the states
    beta    = tau*beta;         % update beta

    betaVect(t) = beta;         % store results over time
    EVect(t)    = -0.5 * sum(sum(W .* (S' * S)));
end

figure(3)
subplot 211
plot(EVect)
title('Energy (mean field annealing)');
xlabel('time');
ylabel('energy');
subplot 212
plot(1./betaVect)
title('beta');
xlabel('time');
ylabel('temperature');




%% Exercise 3

K       = 4;
eta     = 0.3;
tau     = 0.8;
w       = repmat(mean(data,2), 1, K) + randn(2,K) .* repmat(0.1*std(data,0,2),1,K);
t_max   = length(data);
all_w   = cell(t_max);
err     = NaN(1,t_max);

for t = 1:t_max
    
    % annealing
    if t > t_max/4 
        eta = tau*eta;
    end
    
    % compute distances
    dists = NaN(K,length(data));
    for i = 1:K
        dists(i,:) = sum((data - repmat(w(:,i),1,length(data))).^2);
    end
    [mins, idxs] = min(dists);
    random_idx   = randi(length(data));
    w(:,idxs(random_idx)) = w(:,idxs(random_idx)) - eta * (w(:,idxs(random_idx)) - data(:,random_idx));
    all_w{t} = w;
    err(t)  = sum(mins);
end

figure(4)
plots = [1, 10, 20, 30, 50, 500];
for i = 1:length(plots)
    subplot(2,3,i);
    scatter(data(1,:), data(2,:))
    hold on
    plot(all_w{plots(i)}(1,:), all_w{plots(i)}(2,:),'rx')
    hold off
    title(['iteration: ' int2str(plots(i))]);
end

figure(5)
plot(err)