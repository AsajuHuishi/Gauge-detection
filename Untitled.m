img = imread('53_rm.png');
img(:) = 0;
figure;
imshow(img);%全黑
hold on     % 画一条线
x = [100 200];
y = 2*x;
mask = plot(x,y);
imwrite(uint8(mask),'m.png') %保存为‘m.png’