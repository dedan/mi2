 
 close all
 clear
 clc
 
 sz = get(0,'ScreenSize');
 
%% Exercise 1 - Simulated annealing

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
    random_node = randint(1,1,N)+1;     % select random node
    
    % compute local energy of this node
    E = -0.5 *S(random_node) * (W(random_node,:)*S');
    
    % compute probability to change (random_node)
    p = 1 / (1-exp(-beta*2*E));
    if(rand < p)
        S(random_node) = -S(random_node);
    end
    beta        = tau*beta;     % update beta
    
    % store results over time
    betaVect(t) = beta;         
    EVect(t)    = -0.5 * sum(sum(W .* (S' * S)));
end
final_S = S;

fig = figure(1);
set(fig,'PaperPositionMode','auto');
annotation('textbox','String','Simulated annealing','Position',[0.45 1 0.5 0],'LineStyle','none');
subplot 211
plot(EVect)
xlabel('time');
ylabel('energy');
subplot 212
plot(1./betaVect)
xlabel('time');
ylabel('temperature');
print -depsc report/ex1_annealing_et.eps

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

fig = figure(2);
set(fig,'PaperPositionMode','auto');
set(fig,'Position',[0 sz(4) sz(3)*0.5 sz(4)]);
subplot 611;
bar(E_all);
ylabel('E(S)');

for i=1:length(betas)
    subplot(6,1,i+1);
    bar(P_all(i,:));
    ylabel('P(S)');
    title(sprintf('beta = %.1f',betas(i)));
end
xlabel('state');
print -depsc report/ex1_annealing_bar.eps

%% Exercise 2 - Mean-fields annealing

N = 6;                          % number of nodes
beta = 0.1;                     % inverse of temperature
tau  = 1.1;                     % annealing factor
tmax = 100;                     % maximum time

betaVect    = zeros(1,tmax);    % to save betas over time
EVect       = zeros(1,tmax);    % to save total network energy over time
S = S_init;


% iterate over time
for t=1:tmax
    random_node = randint(1,1,N)+1;    % select random node

    e = W(random_node,:) * S';         % compute the mean field
    S(random_node) = tanh(beta * e');  % update the states
    beta = tau*beta;                   % update beta

    % store results over time
    betaVect(t) = beta;         
    EVect(t)    = -0.5 * sum(sum(W .* (S' * S)));
end


fig = figure(3)
set(fig,'PaperPositionMode','auto');
annotation('textbox','String','Mean-field annealing','Position',[0.45 0.99 0.5 0],'LineStyle','none');
subplot 211
plot(EVect)
xlabel('time');
ylabel('energy');
subplot 212
plot(1./betaVect)
xlabel('time');
ylabel('temperature');
print -depsc report/ex2_mean_fields_et.eps



%% Exercise 3 - Online K-means clustering 

load data.dat

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
    random_idx   = randint(1,1,length(data))+1;
    w(:,idxs(random_idx)) = w(:,idxs(random_idx)) - eta * (w(:,idxs(random_idx)) - data(:,random_idx));
    all_w{t} = w;
    err(t)  = sum(mins);
end
err = err / (2*length(data));

fig =figure(4);
set(fig,'PaperPositionMode','auto');
set(fig,'Position',[0 sz(4)*0.5 sz(3)*0.7 sz(4)*0.5]);
plots = [1, 10, 20, 50, 100, 500];
for i = 1:length(plots)
    subplot(2,3,i);
    scatter(data(1,:), data(2,:))
    hold on
    plot(all_w{plots(i)}(1,:), all_w{plots(i)}(2,:),'r.')
    hold off
    title(['iteration: ' int2str(plots(i))]);
end
print -depsc report/ex3_kmeans_steps.eps


fig = figure(5);
set(fig,'PaperPositionMode','auto');
plot(err)
title('K-means error function')
xlabel('iteration')
ylabel('error')
print -depsc report/ex3_kmeans_error.eps

%% Exercise 4 - Soft K-means with fixed betas

load data.dat
K       = 8;
w_init  = repmat(mean(data,2), 1, K) + randn(2,K) .* repmat(0.5*std(data,0,2),1,K);
w_old   = w_init;
w_new   = NaN(size(w_old));
m = NaN(K,length(data));
gamma   = 0.01;
betas   = 0.2:0.2:20;
all_w   = cell(length(betas));

for i = 1:length(betas)
    beta = betas(i);
    stop = false;
    while ~stop
        for q = 1:K;
            m(q,:) = exp(- 0.5*beta*sum((data - repmat(w_old(:,q),1,length(data))).^2));
        end
        m = m ./ repmat(sum(m),K,1);
        
        for q = 1:K;
            w_new(1,q) = m(q,:)*data(1,:)';
            w_new(2,q) = m(q,:)*data(2,:)';
            w_new(:,q) = w_new(:,q) ./ sum(m(q,:));
        end
        all_w{i} = w_new;
        stop = sum(sum((w_new - w_old).^2) < ones(1,K)*gamma) == K;
        w_old = w_new;

    end
end

fig =figure(6);
set(fig,'PaperPositionMode','auto');
set(fig,'Position',[0 sz(4)*0.5 sz(3)*0.5 sz(4)*0.5]);
annotation('textbox','String','Soft K-means (fixed beta)','Position',[0.45 0.99 0.5 0],'LineStyle','none');
scatter(data(1,:), data(2,:))
hold on
plot(w_init(1,:), w_init(2,:),'g.')
for i = 1:length(betas)
    plot(all_w{i}(1,:), all_w{i}(2,:),'r.');
end
legend('data','initial prot.','final prot.')
print -depsc report/ex4_soft_kmeans.eps

fig = figure(7);
set(fig,'PaperPositionMode','auto');
for i = 1:length(betas)
    hold on
    plot(ones(1,K)*betas(i),all_w{i}(1,:),'.');
end
hold off
title('prototypes'' first coordinate as a function of beta')
xlabel('beta')
ylabel('prototypes'' first coordinate')
print -depsc report/ex4_soft_kmeans_protbeta.eps

%% Exercise 4  Soft-K-means with annealing

load data.dat

Ks = [ 4 6 8];
beta_init = 0.2;
beta_max = 100;
tau = 1.1;
all_w_final= cell(length(Ks));
all_w_init = cell(length(Ks)); 

for i = 1:length(Ks)
    
    K = Ks(i);
    w_init  = repmat(mean(data,2), 1, K) + randn(2,K) .* repmat(0.5*std(data,0,2),1,K);
    w_old   = w_init;
    w_new   = NaN(size(w_old));
    m = NaN(K,length(data));

    all_w_init{i} = w_init; 
    beta = beta_init;
    
    while beta <= beta_max
        stop = false;
        while ~stop
            for q = 1:K;
                m(q,:) = exp(- 0.5*beta*sum((data - repmat(w_old(:,q),1,length(data))).^2));
            end
            m = m ./ repmat(sum(m),K,1);

            for q = 1:K;
                w_new(1,q) = m(q,:)*data(1,:)';
                w_new(2,q) = m(q,:)*data(2,:)';
                w_new(:,q) = w_new(:,q) ./ sum(m(q,:));
            end
            stop = sum(sum((w_new - w_old).^2) < ones(1,K)*gamma) == K;
            w_old = w_new;

        end 
        beta = tau*beta;
    end
    all_w_final{i} = w_new;
end

fig = figure(8);
set(fig,'PaperPositionMode','auto');
set(fig,'Position',[0 sz(4)*0.5 sz(3)*0.85 sz(4)*0.4]);
%annotation('textbox','String','Soft K-means (with annealing)','Position',[0.45 0.99 0.5 0],'LineStyle','none');
for i = 1:length(Ks)
    subplot(1,length(Ks),i)
    scatter(data(1,:), data(2,:))
    hold on
    plot(all_w_init{i}(1,:), all_w_init{i}(2,:),'g.');
    plot(all_w_final{i}(1,:), all_w_final{i}(2,:),'r.');
    legend('data','initial prototypes','final prototypes')
    title(sprintf('K = %d',Ks(i)));
end
print -depsc report/ex4_soft_kmeans_annealing.eps
