function patches = get_random_patches(prefix,num_patches,patch_size)

    path = 'imgpca';                                % path to the image folder
    types = ['n','b','d','t'];                      % type of images
    max_imgs = [13, 10 ,12, 14 ];                   % maximum nr. of images for each type
    hps = round(patch_size/2);                      % half of the patch size 
    patches = zeros(patch_size^2,num_patches);      % final matrix of patches
    max_img_type = max_imgs(types==prefix);         % max num of imgs for the requested type
    
    % check if the image type is valid
    if sum(types==prefix) == 1
        
        patch_per_img = floor(num_patches/max_img_type);    % nr. patches per image
    
        % iterate over all images
        for i = 1:max_img_type
            img_file = fullfile(path,strcat(prefix,num2str(i),'.jpg'));  % image file
            img = imread(img_file);                                      % read image
            [width, height] = size(img);                                 % image size
        
            % extract the needed number of patches from this image
            for j = 1:patch_per_img
                center_x = randi(1,1,width-patch_size)+hps;     % x center of the patch
                center_y = randi(1,1,height-patch_size)+hps;    % y center of the patch
	
                patch = img(center_x-hps+1:center_x+hps, center_y-hps+1:center_y+hps);  % extract patch
                patches(:,i) = reshape(patch,patch_size^2,1);                           % reshape and store the patch
            end
        end
    else
        error('image type not valid')
    end
end