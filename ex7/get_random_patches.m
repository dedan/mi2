function patches = get_random_patches(prefix,num_patches,patch_size)

    if exist('randi')
        randintfun = @randi;
    else
        randintfun = @randint;
    end
    
    path = 'images';                                % path to the image folder
    types = ['n','b','d','t'];                      % type of images
    max_imgs = [13, 10 ,12, 14 ];                   % maximum nr. of images for each type
    hps = round(patch_size/2);                      % half of the patch size 
    patches = zeros(patch_size^2,num_patches);      % final matrix of patches
    max_img_type = max_imgs(types==prefix);         % max num of imgs for the requested type
    
    % check if the image type is valid
    if sum(types==prefix) == 1
        
        patch_per_img = ceil(num_patches/max_img_type);    % nr. patches per image
        count = 1;
        
        % iterate over all images
        for i = 1:max_img_type
            img_file = fullfile(path,strcat(prefix,num2str(i),'.jpg'));  % image file
            img = imread(img_file);                                      % read image
            [width, height] = size(img);                                 % image size
        
            % extract the needed number of patches from this image
            for j = 1:patch_per_img
                center_x = randintfun(1,1,width-patch_size)+hps;     % x center of the patch
                center_y = randintfun(1,1,height-patch_size)+hps;    % y center of the patch
	
                patch = img(center_x-hps+1:center_x+hps, center_y-hps+1:center_y+hps);  % extract patch
                patches(:,count) = reshape(patch,patch_size^2,1);                       % reshape and store the patch
                if count == num_patches
                    break;
                end
                count = count+1;
            end
        end
    else
        error('image type not valid')
    end
end