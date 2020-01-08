% 'edge_extraction.m'+'line_fitting.m'+'elipse_fitting.m'
pic = 6;
path = strcat('orig\orig\',num2str(pic),'.jpg');
img = imread(path);
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

imshow([emx, emy, em], [])
savepath = strcat('orig-operation\',num2str(pic),'.png');
imwrite(em, savepath);

path = strcat('orig-operation\',num2str(pic),'.png');% em
%  X = double(imread('33.png'));% Ô­À´µÄ
X = 255 * double(imread(path));

[M, N] =  size(X);
indx = 1:N;
indy = 1:M;
posidx_curvpt = find(X==255);
n_curvpt = length(posidx_curvpt);

%line fitting
n_ptpair = round(0.03*n_curvpt);
rng(10);
idx_ptpair = [randi(n_curvpt, n_ptpair, 1),randi(n_curvpt, n_ptpair, 1)]; 
X_disp = zeros(M, N, 3);
X_disp(:,:, 1) = X;
X_disp(:,:, 2) = X;
X_disp(:,:, 3) = X;
X_disp(M*N + posidx_curvpt(idx_ptpair(:,1))) = 0;
X_disp(2*M*N + posidx_curvpt(idx_ptpair(:,1))) = 0;
X_disp(M*N + posidx_curvpt(idx_ptpair(:,2))) = 0;
X_disp(2*M*N + posidx_curvpt(idx_ptpair(:,2))) = 0;

[pt1y, pt1x] = ind2sub([M, N], posidx_curvpt(idx_ptpair(:,1)));
[pt2y, pt2x] = ind2sub([M, N], posidx_curvpt(idx_ptpair(:,2)));
[ptally, ptallx] = ind2sub([M, N], posidx_curvpt);

k = (pt1y - pt2y)./(pt1x - pt2x + eps);
b = pt2y - k.*pt2x;
  
figure(100)
imshow(uint8(X_disp));
for i=1:n_ptpair
    line(indx, k(i).*indx + b(i),  'Color', 'g')
end

% mask_verline = k>M;
distmap = zeros(n_curvpt, n_ptpair);
for i=1:n_ptpair
    %sampling
    distmap(:, i) = abs( ptally - k(i)*ptallx -b(i) )/sqrt(k(i)^2 + 1);
end
simple_voting = 1;
if simple_voting
    %% simple voting
    distvote = sum(distmap <3, 1);
    maxvote = max(distvote(:));
    is_voted = distvote(:) > 0.4*maxvote;
else
%% voting with line merging
    maskvote = double(distmap<3);
    corr_lines = maskvote'*maskvote;
    autocorr_lines = diag(corr_lines);
    figure(100)
    imshow(uint8(X_disp));
    for i=1:n_ptpair
        line(indx, k(i).*indx + b(i),  'Color', 'g')
    end
    colorpool = ['r', 'g', 'w'];
    for j = 1:10
        [max_fit, ind_max_fit] = max(autocorr_lines);
       
        corr_w_max = corr_lines(ind_max_fit, :);
        ind_line_cluster = find(autocorr_lines>0.9*max_fit);
         disp([max_fit, length(ind_line_cluster)])
        for i=1:length(ind_line_cluster)
            line(indx, k(ind_line_cluster(i)).*indx + b(ind_line_cluster(i)),  'Color', colorpool(mod(j, 3)+1)  );
        end
        autocorr_lines(ind_line_cluster) = 0;
        pause
    end
    
end

X_rm = X;
figure(101)
imshow(uint8(X_disp));
for i=1:n_ptpair
    if is_voted(i)
        line(indx, k(i).*indx + b(i),  'Color', 'g');
         X_rm(ind2sub([M, N], posidx_curvpt(distmap(:, i)<5))) = 0;
    end
end
figure(102)
imshow(uint8(X_rm));
savepath = strcat('remove-line\',num2str(pic),'.png');
imwrite(X_rm, savepath);

if(0)
%circle fitting
X = X_rm;
posidx_curvpt = find(X==255);
n_curvpt = length(posidx_curvpt);
n_pttriple = round(0.01*n_curvpt);
rng(10);
idx_pttripler = [randi(n_curvpt, n_pttriple, 1),randi(n_curvpt, n_pttriple, 1),randi(n_curvpt, n_pttriple, 1)]; %#ok<DRNDINT>
[pt1y, pt1x] = ind2sub([M, N], posidx_curvpt(idx_pttripler(:,1)));
[pt2y, pt2x] = ind2sub([M, N], posidx_curvpt(idx_pttripler(:,2)));
[pt3y, pt3x] = ind2sub([M, N], posidx_curvpt(idx_pttripler(:,3)));
[ptally, ptallx] = ind2sub([M, N], posidx_curvpt);
y1 = pt1y; y2 = pt2y; y3 = pt3y;
x1 = pt1x; x2 = pt2x; x3 = pt3x; 
x2 = x2 - x1; x3 = x3 - x1;
y2 = y2 - y1; y3 = y3 - y1;
m = 2*(x2.*y3 - y2.*x3);
cx = (x2.^2 .* y3 - x3.^2 .*y2 + y2.*y3 .* (y2-y3))./m;
cy = (x2 .* x3 .* (x3-x2) - y2.^2 .* x3 + x2 .* y3.^2)./m;
r = sqrt(cx.^2 + cy.^2);
cx = cx + x1;
cy = cy + y1;


X_disp = zeros(M, N, 3);
X_disp(:,:, 1) = X;
X_disp(:,:, 2) = X;
X_disp(:,:, 3) = X;
X_disp(M*N + posidx_curvpt(idx_pttripler(:,1))) = 0;
X_disp(2*M*N + posidx_curvpt(idx_pttripler(:,1))) = 0;
X_disp(M*N + posidx_curvpt(idx_pttripler(:,2))) = 0;
X_disp(2*M*N + posidx_curvpt(idx_pttripler(:,2))) = 0;
X_disp(M*N + posidx_curvpt(idx_pttripler(:,3))) = 0;
X_disp(2*M*N + posidx_curvpt(idx_pttripler(:,3))) = 0;
ph = figure(103)
imshow(uint8(X_disp));
for i=1:n_pttriple
    circle(ph, cx(i), cy(i), r(i))
end

distmap = zeros(n_curvpt, n_pttriple);
for i=1:n_pttriple
    %sampling
    distmap(:, i) = abs( sqrt( (ptallx - cx(i)).^2  + (ptally - cy(i)).^2) - r(i));
end
distvote = sum(distmap <3, 1);
maxvote = max(distvote(:));
is_voted = distvote(:) > 0.8*maxvote;
is_voted = is_voted.* (r<min(M, N)/2);


ph = figure(104)
imshow(uint8(X_disp));
for i=1:n_pttriple
    if is_voted(i)
        circle(ph, cx(i), cy(i), r(i))
    end
end

end

% 'ellipse_fitting.m' Author:Zeng Xinyang,TJU ,19-04-18
% need function with 'fit_ellipse_state.m','fit_ellipse_state.m','fit_ellipse_state2.m'
% TODO: detect ellipse in image, which is marked in yellow curve.
% Notice: this code is based on Probability sampling, the results may vary from run to run. 
% Please run it again if not ideal in one time.
savepath = strcat('./result425/');
mkdir(savepath)

path = strcat('remove-line\',num2str(pic),'.png')
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
    to = to +1;    
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
% figure();
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

close all
clearvars -EXCEPT Indicator indi




