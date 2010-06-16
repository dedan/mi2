
close all
clear

%% part 1 init

% load the sounds
load sounds/sound1.dat
load sounds/sound2.dat
s = [sound1'; sound2'];

% settings 
N = 2;

% create mixing matrix and mix sources
A = rand(N);
x = A*s;

% permute columns of data matrix
x_perm = x(:,randperm(length(x)));

% correlation coefficients
for i = 1:N
    c = corrcoef(s(i,:)', x(i,:)');
    disp(['correlation signal ' int2str(i) ': ' num2str(c(2,1))]);
end

% center the data
x = bsxfun(@minus, x, mean(x,2));

% initialize unmixing matrix w randomly
w = rand(N);


%% part 2 optimization



% soundsc(sound1)
% soundsc(sound2)
