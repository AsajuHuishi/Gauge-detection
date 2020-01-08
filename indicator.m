%% get precision, recall, f-score from 'ellipse_fitting.m'
% 19.4.23-24
%     imshow(gt_path);
clear
Indicator = [];
for pic = 74:118
    gt_path = strcat('real-bin\',num2str(pic),'.png');
    img_path = strcat('result423\',num2str(pic),'_mask.png');
  
    gt = imread(gt_path);
    img = imread(img_path);
%     figure;
%     imshow(gt)
%     figure;
%     imshow(img_path);
    %% img to red
    image = img(:,:,1);
    [r,c] = size(image);
    image = ones(r,c,3);
    for i = 1:r
        for j = 1:c
            if img(i,j) == 255
                image(i,j,1) = 255;
                image(i,j,2) = 0;
                image(i,j,3) = 0;
            end
        end
    end
%     figure;
%     imshow(image);
    %% gt to green
    gt2 = gt(:,:,1);
    [rr,cc] = size(gt2);
    gt2 = ones(rr,cc,3);
    for ii = 1:rr
        for jj = 1:cc
            if gt(ii,jj) == 255
                gt2(ii,jj,1) = 0;
                gt2(ii,jj,2) = 255;
                gt2(ii,jj,3) = 0;
            end
        end
    end
%     figure;
%     imshow(gt2);
   %% add gt2 and image 
    image = uint8(image);
    gt2 = uint8(gt2);
    whole = imadd(image,gt2);
%     figure;
%     imshow(whole);    
    %% calculate indicators
    gt = gt(:,:,1)>128;
%     figure;
%     imshow(gt);
     img = img(:,:,1)>128;
%     figure;
%     imshow(img);
    gt = gt(:);
    img = img(:);
    a = sum(gt);
    b = sum(img);
    c = sum(gt&img);
    pre = c/b;
    rec = c/a;
    F1 = 2*pre*rec/(pre + rec);
    disp(['image  ',num2str(pic),':',num2str(pre),'    ',num2str(rec),'    ',num2str(F1)]);
    indi = [pic pre rec F1];
    Indicator = [Indicator; indi];
    clearvars -EXCEPT Indicator indi
end
