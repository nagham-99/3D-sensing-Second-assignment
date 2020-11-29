
img1 = double(imread('einstein.jpg'))/255;
img2 = double(imread('Lenna.png'))/255;


img1 = img1+0.03*randn(size(img1));
img2 = img2+0.03*randn(size(img2));
img1(img1<0) = 0; img1(img1>1) = 1;
img2(img2<0) = 0; img2(img2>1) = 1;


w     = 2;       
sigma = [3 0.1]; 


bflt_img1 = bfilter2(img1,w,sigma);
bflt_img2 = bfilter2(img2,w,sigma);


figure(1); clf;
set(gcf,'Name','Grayscale Bilateral Filtering Results');
subplot(1,2,1); imagesc(img1);
axis image; colormap gray;
title('Input Image');
subplot(1,2,2); imagesc(bflt_img1);
axis image; colormap gray;
title('Result of Bilateral Filtering');


figure(2); clf;
set(gcf,'Name','Color Bilateral Filtering Results');
subplot(1,2,1); imagesc(img2);
axis image; colormap gray;
title('Input Image');
subplot(1,2,2); imagesc(bflt_img2);
axis image; title('Result of Bilateral Filtering');
drawnow;

