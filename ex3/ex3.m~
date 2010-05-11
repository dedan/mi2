
im          = imread('testimg.jpg');
sigmas      = [0.05 0.1] * double(max(max(im)));
sample_size = [100 500];

noise   = randn(size(im)) * sigma(1);
imshow(double(im) + noise);

perms = randperm(size(im(:)));
random_sample = im(perms(1:sample_size(1)));

for h = 1:10:50
   kernel = ones(1,h) ./ h;
   c = conv(random_sample, kernel);
   c = c(length(kernel)/2:length(c)-length(kernel)/2);
    
end