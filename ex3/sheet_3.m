clear all
close all
clc

img = imread('testimg.jpg');    % read image
img = im2double(img);           % convert to double and normalize

noiseSTD = [0.05 0.1];          % std for white noise   
sampleSize = [100 500];

imgN = imnoise(img,'gaussian', 0, noiseSTD(1)); % add noise

imgNvec = reshape(imgN, 1, []); % flaten image

sample = randsample(imgNvec, sampleSize(1)); % get random sample set from image

cutImg = imgNvec;

for i = 1:length(sample)
    toBeRemoved = find(cutImg == sample(i));
    cutImg(toBeRemoved(1)) = []; % remove samples used for estimation
end

sigma = 0.1; % sigma for kernel
range = -1:0.01:1;
kernel = [];

for x = range
    kernel = cat(1,kernel,gauss(x,sigma)); % compute kernel
end

figure
plot(range,kernel) % plot kernel

binWidth = 0.01; % width for histogram
binVec = min(sample):binWidth:max(sample);

h = histc(sample,binVec); % compute histogram
h = h./binWidth;
figure
bar(binVec, h, 'histc')

figure
f = conv(kernel,sample); % convolve kernel with sample set
plot(f)

figure
f = conv(kernel,h); % convole kernel with histogram
plot(f)

figure()
f = conv(kernel,cutImg); % convole kernel with img
plot(f)

%imshow(imgN)


