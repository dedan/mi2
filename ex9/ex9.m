%% Exercise 1

close all
clear
clc

N = 6;

% randomly initialize initial states
S = ones(1,N);
S(rand(1,N) > 0.5) = -1;

% randomly initialize weights
W = rand(N,N);
W = W+W';
W = W-diag(diag(W)); 

beta = 0.1;
tau = 1.1;
tmax = 10;

betaVect = zeros(1,tmax);
EVect = zeros(N,tmax);
%Eold = Inf;
for t=1:tmax
    i = randint(1,1,N)+1;
    E = -0.5 *S(i) * (W(i,:)*S');
    p = 1 / (1-exp(-beta*2*E));
    if(rand < p)
        S(i) = -S(i);
    end
    beta = tau*beta;
end