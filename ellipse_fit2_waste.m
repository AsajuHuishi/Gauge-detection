load('sample_mat','sample')
N = 1000;
sample_x = sample(:,1);
sample_y = sample(:,2);
% hold on
% plot(sample_y,sample_x,'o')
% set(gca,'XAxisLocation','top'); 
% set(gca,'YDir','reverse');  
iter = 2000;
count_elli = 0;
distance = zeros(N,iter);
res = zeros(N,iter);
Nk = zeros(11,iter);
fixed = zeros(1,iter);
x_state = zeros(5,1);
y_state = zeros(5,1);
for i = 1:iter
    sam_index = randperm(N,5);
    x = sample_x(sam_index);
    y = sample_y(sam_index);

%     hold on;
%     plot(y,x,'*');
%     set(gca,'XAxisLocation','top'); 
%     set(gca,'YDir','reverse');  
    x2=x.^2;
    y2=y.^2;
    xy=x.*y;
    A=[x2,xy,y2,x,y];
    c = inv(A)*ones(5,1);
    nabla = det([c(1) c(2) c(4);c(2) c(3) c(5);c(4) c(5) -1]);
    delta = det([c(1) c(2); c(2) c(3)]);
    S = c(1) + c(3);
    if(delta>0&&nabla~=0&&nabla*S<0)
        count_elli = count_elli + 1;
        hold on
        [focus1, focus2, long_axis] = fit_ellipse(x,y,1);
        %% check the number sample points in the ellipse
        for row = 1:1000
            distance(row,count_elli) = norm(sample(row,:) - focus1) + norm(sample(row,:) - focus2);
        end
        res(:,count_elli) = distance(:,count_elli) - long_axis;
        Nk(1,count_elli) = sum(abs(res(:,count_elli))<=0.1); % tau
        Nk(2:6,count_elli) = x;
        Nk(7:11,count_elli) = y;    
    end    
end
%% draw 5 ellipse with biggest Nk
img = imread('33_rm.png');
x_n = im2double(img);
figure();
imshow(img)
hold on

Nk(:,all(Nk==0,1)) = [];
Nk1 = sortrows(Nk',1);
for a = 0:4
% a=0
    x_state = Nk1(end-a,2:6)';
    y_state = Nk1(end-a,7:11)';
    fit_ellipse_state(x_state,y_state,1);
end
% Author:Zeng Xinyang,19-04-18
% figure();
% imshow(img)