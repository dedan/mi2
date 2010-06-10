function ex5_gr1

% MI ex?? plotting Laplacian distribution

% Covariance Matrix

C           = [4,-1;-1,2];
m           = [2;-1];
n_samples   = 500;

%Diagonalize C Matrix
[A,B] = eig(C);

L = [laplace_dist(n_samples); laplace_dist(n_samples)];     % create 2d dist
L = sqrt(B) * L;                                            % scale the variances
L = A'*L + repmat(m,[1,length(L)]);                         % introduce correlation

disp('Covariance of sample points:') 
disp(cov(L'))

disp('Mean of sample points:')
disp(mean(L,2)); 

figure(1);
plot(L(1,:),L(2,:),'.');


function L = laplace_dist(N)

%Required formula is:-
% X = mu - b*sign(Y)*log(1-2*abs(Y))

mean_of_uniform = 0.5;
b               = 1/sqrt(2);
uniform         = rand(1,N) - mean_of_uniform;

L = b*sign(uniform).*log(1-2*abs(uniform));



