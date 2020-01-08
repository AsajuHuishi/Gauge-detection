pic = 83;
path = strcat('orig\orig\',num2str(pic),'.jpg');
img = imread(path);
figure;
imshow(img)
img = double(rgb2gray(img));
img2 = medfilt2(img, [21 21]);
hx = [-1 0 1; -2 0 2; -1 0 1];
hy = [-1 -2 -1; 0 0 0; 1 2 1];

imgx = imfilter(img2, hx, 'same');
imgy = imfilter(img2, hy, 'same');
mask = ones(size(img));
Tx = 150;
Ty = 100;
mask(1:Ty, :) = 0;
mask(end-Ty:end, :) = 0;
mask(:,1:Tx) = 0;
mask(:, end-Tx:end) = 0;
imgx = mask.*imgx;
imgy = mask.*imgy;
imgxy = sqrt(imgx.^2 + imgy.^2)/sqrt(2);
tau = 255/5;
emx = abs(imgx) > tau;
emy = abs(imgy) > tau;
em =  imgxy > tau;
% imshow([emx, emy, em], [])
savepath = strcat('orig-operation\',num2str(pic),'.png');
imwrite(em, savepath);

B=[0 0 1 0 0; % ≈Ú’Õ
    0 1 1 1 0;
    1 1 1 1 1;
    0 1 1 1 0;
    0 0 1 0 0];
em = imdilate(em, B);
figure;
imshow(em);
disppath = strcat('orig-operation2\',num2str(pic),'.png');
imwrite(em, disppath);