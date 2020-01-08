% ͼ1���Ա�����������Բ����Բ���
% Zxy, 19-5-5��v2 5-8��
%% ��ȡ ���������� mask
pic = 42;% 74 42 14
houghpath = strcat('hough\hough\',num2str(pic),'.png');
path1 = strcat('template\template\ccoeff\',num2str(pic),'.png');
path2 = strcat('template\template\ccoeff_normed\',num2str(pic),'.png');
path3 = strcat('template\template\ccorr\',num2str(pic),'.png');
path4 = strcat('template\template\ccorr_normed\',num2str(pic),'.png');
path5 = strcat('template\template\sqdiff\',num2str(pic),'.png');
path6 = strcat('template\template\sqdiff_normed\',num2str(pic),'.png');

im = imread(houghpath);    
im1 = imread(path1);
im2 = imread(path2);
im3 = imread(path3);
im4 = imread(path4);
im5 = imread(path5);
im6 = imread(path6);

im = im2double(im);
im1 = im2double(im1);           % coeff
im2 = im2double(im2);           % coeff norm
im3 = im2double(im3);           % ccorr
im4 = im2double(im4);           % ccorr norm
im5 = im2double(im5);           % sqdiff
im6 = im2double(im6);           % sqdiff norm

% 6����1
im = im4         ;

% ת��ֵͼ��  
bw = im2bw( im );  
  % �߽���  
contour = edge(bw ,'canny'); 
%% ��ȡ our
ourpath = strcat('./Ŀǰ��result/',num2str(pic),'_mask.png');
im2 = imread(ourpath);
im2 = im2double(im2);

% ת��ֵͼ��  
bw2 = im2bw( im2 );                         
% �߽���  
contour2 = edge(bw2 ,'canny');  
s = zeros(1080,1920);
for i = 1:1080
    for j = 1:1920
        if contour2(i,j) == 1
            s(i,j,1) = 1;
            s(i,j,2) = 1;
            s(i,j,3) = 0;
        end
    end
end
B=[0 0 1 0 0; % ����
    0 1 1 1 0;
    1 1 1 1 1;
    0 1 1 1 0;
    0 0 1 0 0];
s = imdilate(s, B);
figure;
imshow(s)
disppath1 = strcat('contour2\',num2str(pic),'.png');
imwrite(s, disppath1);
%% hough��ԭͼ
path = strcat('orig\orig\',num2str(pic),'.jpg');
orig = imread(path);
orig = im2double(orig);
figure;
imshow(orig);
origx = orig;
origy = orig;

for i = 1:1080
    for j = 1:1920
        if contour(i,j) == 1
            origx(i,j,1) = 1;
            origx(i,j,2) = 1;
            origx(i,j,3) = 0;
        end
    end
end
origx = imdilate(origx, B);
figure;
imshow(origx)

%% our��ԭͼ
for i = 1:1080
    for j = 1:1920
        if contour2(i,j) == 1
            origy(i,j,1) = 1;
            origy(i,j,2) = 1;
            origy(i,j,3) = 0;
        end
    end
end
B2=[0 0 1 0 0; % ����
    0 1 1 1 0;
    1 1 1 1 1;
    0 1 1 1 0;
    0 0 1 0 0];
origy = imdilate(origy, B2);
figure;
imshow(origy)
            
disppath2 = strcat('final2\',num2str(pic),'.png');
imwrite(origy, disppath2);