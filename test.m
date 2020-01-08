clear
img = imread('71_rm.png');
  x_n = im2double(img);
[pos_x, pos_y] = find(x_n==1);
pos = [pos_x pos_y];

N = size(pos,1);
Np = 6;
sam_index = randperm(N,Np);
sample = pos;
sample_x = sample(:,1);
sample_y = sample(:,2);
x = sample_x(sam_index);
y = sample_y(sam_index);

figure
imshow(img)
hold on
[f1, f2, lang] = fit_ellipse( x,y,1 );

