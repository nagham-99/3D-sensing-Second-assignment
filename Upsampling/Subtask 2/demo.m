
addpath('Inputs');

  scale = 1; 
                
guide_s = double(rgb2gray(imread('art_color.png')));
guide_s = guide_s/max(guide_s(:));
depth_gt = double(imread('art.png'));

% make the input LR depth map.
h = fspecial('gaussian', 2^scale, 2^scale);
tmp = imfilter(depth_gt, h, 'replicate');
depth = tmp(1:2^scale:end, 1:2^scale:end);
d_max = max(depth(:));
depth = depth/d_max;

if ~exist('Results', 'dir')
    mkdir('Results');
end
savepath = sprintf('Results\\result.png', ...
                    size(depth), size(depth_gt));

% perform x2 upsampling operation iteratively.
for i = 1:scale
    h = fspecial('gaussian', 2^(scale-i), 2^(scale-i));
    guide_tmp = imfilter(guide_s, h, 'replicate');
    guide = guide_tmp(1:2^(scale-i):end, 1:2^(scale-i):end);
    depth = resint(depth, guide);
end

imshow(uint8(depth*d_max))
imwrite(uint8(depth*d_max), savepath);