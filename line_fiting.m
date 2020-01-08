pic = 83;
path = strcat('orig-operation\',num2str(pic),'.png');% em
%  X = double(imread('33.png'));% ‘≠¿¥µƒ
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
    distvote = sum(distmap <4, 1);
    maxvote = max(distvote(:));
    is_voted = distvote(:) > 0.5*maxvote;
else
%% voting with line merging
    maskvote = double(distmap<3);
    corr_lines = maskvote'*maskvote;
    autocorr_lines = diag(corr_lines);
    figure(100)
    imshow(uint8(X_disp));
    for i=1:n_ptpair
         line(indx, k(i).*indx + b(i),  'Color', 'g')
%     X_disp(indx, k(i).*indx + b(i))=
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
B=[0 0 1 0 0; % ≈Ú’Õ
    0 1 1 1 0;
    1 1 1 1 1;
    0 1 1 1 0;
    0 0 1 0 0];
X_disp = imdilate(X_disp, B);
figure(101)
imshow(uint8(X_disp));

for i=1:n_ptpair
    if is_voted(i)
        line(indx, k(i).*indx + b(i),  'Color', 'g');
         X_rm(ind2sub([M, N], posidx_curvpt(distmap(:, i)<5))) = 0;
    end
end
savepath = strcat('remove-line\',num2str(pic),'.png');
imwrite(X_rm, savepath);

B=[0 0 1 0 0; % ≈Ú’Õ
    0 1 1 1 0;
    1 1 1 1 1;
    0 1 1 1 0;
    0 0 1 0 0];
X_rm = imdilate(X_rm, B);
figure(102)
imshow(uint8(X_rm));
disppath = strcat('remove-line2\',num2str(pic),'.png');
imwrite(X_rm, disppath);
%% waste
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




