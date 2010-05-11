%%% Bayesian linear regression example: polynomial curve fitting
clear all
close all

% define target function
h = inline('-0.25*x.^2 + 0.75*x.^3','x');
xrange = [-1 1];
beta = (1/0.2)^2;   % accuracy of the target function (noise level)

% store target weights for visualization
tw = [-0.25 0.75];

% width of prior distribution p(w)
ALPHA = [1e-2 1 1e2];

% training sample
nSamples = 50;
x = rand(nSamples,1)*diff(xrange) + xrange(1);   
% calculate noise: beta^-1 = sigma^2
noise = randn(nSamples,1)*1/sqrt(beta);   
t = h(x) + noise;    % noisy target 

% define basis functions
N = 2; % number of base functions
for n = 1:N
    f{n} = inline(sprintf('x.^%i',n+1),'x'); % {x^2, x^3}
end

% define raster for x within the interval xrange
xt = linspace(xrange(1),xrange(2),101);
phit = calc_featurevector(f,xt);   % for plotting fits
% define raster for w for plotting later
wrange = [-1 1];
w1g = linspace(wrange(1),wrange(2),41);
w2g = linspace(wrange(1),wrange(2),41);
[W1,W2] = meshgrid(w1g,w2g);
w = [W1(:) W2(:)]';

% loop over different alphas
for nalpha = 1:length(ALPHA)
    alpha = ALPHA(nalpha);
    

    % intialize the prior distribution p(w)
    m0 = zeros(N,1);
    S0_inv = eye(2) / alpha;   %EXERCISE
    S0 = inv(S0_inv);
      
    c = 1;  X = [];  T = [];
    % data points from sample set arriving sequentially 
    for i = 1:nSamples           
        % calculate values of the basis functions
        phi(i,:) = calc_featurevector(f,x(i));

        % update posterior p(w|t)
        % HOMEWORK (Bishop Ch. 3.3.1)
        % calculate S_inv and m using beta, phi, S0_inv, m0, t
%        S_inv = ......;   %EXERCISE
        R = inv(chol(S_inv));  % inversion of ...
        S = R*R';              % matrix S_inv: S = S_inv^-1
%        m = ......;   %EXERCISE

        % update prior distribution (=current posterior distribution)
        S0_inv = S_inv;
        S0 = S;
        m0 = m;
            
        % plotting for certain sample sizes        
        if ~isempty(find([1 5 25]==i, 1))
            figure(nalpha)
            subplot(3,3,c+1)
            plot(x(1:i),t(1:i),'bo')  % plot samples

            % plot 5 curves using w sampled from the current 
            % posterior distribution p(w|t)
            for j = 1:5
%                w_sample = ......;   %EXERCISE (e.g. mvnrnd)
%                y = ......;   %EXERCISE
                hold on
                plot(xt,y,'r')
            end
            xlim(xrange)
            ylim([-1 1])
                
            % plot posterior p(w|t)
            subplot(3,3,c)
            [posterior] = calc_gauss(w(1:2,:),m(1:2),S(1:2,1:2));
            imagesc(w1g,w2g,posterior)
            hold on
            plot(tw(1),tw(2),'w+')
            title(['posterior ' int2str(i)])
            axis tight
            c = c+3;
        end
    end
          
    % plot target function + fitting curve using the MAP weight vector
%    y = ......;   %EXERCISE

    % plot prior p(w)
    subplot(3,3,3)
    [prior] = calc_gauss(w(1:2,:),zeros(N,1),alpha*eye(N));
    imagesc(w1g,w2g,prior)
    hold on
    plot(tw(1),tw(2),'w+')
    title('initial prior')
    axis tight
    subplot(3,3,[6 9])
    hold on
    plot(x,t, 'bo')
    plot(xt,h(xt), 'k')
    plot(xt,y, 'm')
    legend(['N_{data}=' int2str(nSamples)],'t_{noiseless}','w_{MAP} fit',2)
    xlim(xrange); ylim([-1 1.5]);
end
