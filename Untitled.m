img = imread('53_rm.png');
img(:) = 0;
figure;
imshow(img);%ȫ��
hold on     % ��һ����
x = [100 200];
y = 2*x;
mask = plot(x,y);
imwrite(uint8(mask),'m.png') %����Ϊ��m.png��