% 'ellipse_fitting.m' Author:Zeng Xinyang,TJU ,19-04-18
% need function with 'fit_ellipse_state.m','fit_ellipse_state.m','fit_ellipse_state2.m'
% TODO: detect ellipse in image, which is marked in yellow curve.
% Notice: this code is based on Probability sampling, the results may vary from run to run. 
% Please run it again if not ideal in one time.
clear all
savepath = strcat('./result425/');
mkdir(savepath)
for pic = 83
    path = strcat('remove-line\',num2str(pic),'.png')
%     path2 = strcat('orig-operation\',num2str(pic),'.png');
    img = imread(path);
    x_n = im2double(img);
    figure;
    imshow(img)
    %% parameters
    percent = 1;
    tau = 0.2;
    %% find white pixels and get N of them as samples
    [pos_x, pos_y] = find(x_n==1);
    pos = [pos_x pos_y];
    index = 1:size(pos,1);
    [img_masuku, num_masuku] = bwlabel(img,8);% get connected set

    % hold on
    % plot(pos_y,pos_x,'*')
    % set(gca,'XAxisLocation','top'); 
    % set(gca,'YDir','reverse');  
    % cindex = randperm(numel(index));
    % get samples randomly from index
    N = size(pos,1);
    iter = floor(percent*N);
    Np = 5;% points to get ellipse
    % b = index(cindex(1:N));
    sample = pos;% N samples
    %% find Np(>=5) points from samples to get an ellipse
    sample_x = sample(:,1);
    sample_y = sample(:,2);
    % load('sample_mat','sample')
    % hold on
    % plot(sample_y,sample_x,'o')
    % set(gca,'XAxisLocation','top'); 
    % set(gca,'YDir','reverse');

    % initiate
    count_elli = 0;
    distance = zeros(N,iter);
    res = zeros(N,iter);
    Nk = zeros(2*Np+1,iter);
    x_state = zeros(Np,1);
    y_state = zeros(Np,1);
    for i = 1:iter
        sam_index = randperm(N,Np);
        x = sample_x(sam_index);
        y = sample_y(sam_index);
        % check whether sample points are in same connected set
%         Index_masuku = [];
% 
%         index_masuku = diag(img_masuku(x,y));
%         Index_masuku = [Index_masuku index_masuku];
% 
%         same = length(Index_masuku) - length(unique(Index_masuku));
%         if(same~=0)
%             continue
        
        % draw sets of ellipses
        hold on
        [focus1, focus2, long_axis] = fit_ellipse(x,y,1);
        % check whether their focuses are in the image
        if focus1(1)>1080||focus1(2)>1920||focus2(1)>1080||focus2(2)>1920
            continue;
        elseif focus1(1)<0||focus1(2)<0||focus2(1)<0||focus2(2)<0
            continue
        end
        %% get the number of inner points(Nk)
        if focus1 == [0 0];
            continue
        end
        count_elli = count_elli + 1;
        for row = 1:N
            distance(row,count_elli) = norm(sample(row,:) - focus1) + norm(sample(row,:) - focus2);
        end
        res(:,count_elli) = distance(:,count_elli) - long_axis;
        Nk(1,count_elli) = sum(abs(res(:,count_elli))<=tau); % tau

        Nk(2:1+Np,count_elli) = x;
        Nk(2+Np:1+2*Np,count_elli) = y;    
    end    
    %% draw ellipse with biggest 5 Nk
    % vote
    figure();
    imshow(img)
    hold on

    Nk(:,all(Nk==0,1)) = [];
    Nk1 = sortrows(Nk',1);

    Content = [];
    to = 0;

    for a = 0:4
        x_state = Nk1(end-a,2:1+Np)';
        y_state = Nk1(end-a,2+Np:1+2*Np)';
        [focus11, focus22, lang_axis] = fit_ellipse_state(x_state,y_state,1);
        % check their focuses are in the image
        if focus11(1)>1080||focus11(2)>1920||focus22(1)>1080||focus22(2)>1920
            continue;
        elseif focus11(1)<0||focus11(2)<0||focus22(1)<0||focus22(2)<0
            continue
        end
        to = to + 1;    
        for row = 1:N
            dis2(row,a+1) = norm(sample(row,:) - focus11) + norm(sample(row,:) - focus22);
        end
        residue(:,a+1) = dis2(:,a+1) - lang_axis;
        dex = find(abs(residue(:,a+1))<=tau);
        chose = sample(dex,:);
        Content = [Content;chose];
    %     if(size(chose,1)>=5) % only >=5 can get ellipse,otherwise formid
    %         x_state2 = chose(:,1);
    %         y_state2 = chose(:,2);
    %     end
   
end
    
    
    %% 
     figure();
    % 
    x_state2 = Content(:,1);
    y_state2 = Content(:,2);

    hold on
    [f1,f2,longaxis] = fit_ellipse_state2(x_state2,y_state2,1);
    savepath = strcat('./result425/',num2str(pic),'_mask.png');
%%     get outside black point
    figure;
    mask = zeros(1080,1920);

    for xx = 1:1080
        for yy = 1:1920
            p = [xx yy];
            entfernt = norm(p-f1,2)+norm(p-f2,2);
            if entfernt <= longaxis
                mask(xx,yy) = 1;
            end
        end
    end
    imshow(mask);
    imwrite(mask,savepath);

    % check whether their focuses are in the image
    if f1(1)>1080||f1(2)>1920||f2(1)>1080||f2(2)>1920
        warning('ellipse_fitting: ellipse detected might be wrong. Please run the code again or check the image');
    elseif f1(1)<0||f1(2)<0||f2(1)<0||f2(2)<0
        warning('ellipse_fitting: ellipse detected might be wrong. Please run the code again or check the image');
    end
    %% get precision, recall, f-score 
    Indicator = [];

    gt_path = strcat('real-bin\',num2str(pic),'.png');
    img_path = strcat('result425\',num2str(pic),'_mask.png');
  
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

%     clearvars -EXCEPT Indicator indi
end
