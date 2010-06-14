
close all
clear


load sounds/sound1.dat
load sounds/sound2.dat
s = [sound1'; sound2'];

% soundsc(sound1)
% soundsc(sound2)

N = 2;

A = rand(N);
x = A*s;
x = x(:,randperm(length(x)));