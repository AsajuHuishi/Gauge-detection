% 
% x=[4.5596 5.0816 5.5546 5.9636 6.2756]';
% y=[0.8145 1.3685 1.9895 2.6925 3.5265]';
% figure();
% [m,n,a] = fit_ellipse(x,y,1);
% hold on
% % plot(m(2),m(1),'*')
% % plot(n(2),n(1),'*')
% % set(gca,'XAxisLocation','top'); 
% % set(gca,'YDir','reverse'); 
% pos_x = [2 3 5 6]';
% pos_y = [21 32 7 8]';
% pos = [pos_x pos_y]
% % focus1 = [22 33];
% % focus2 = [44 11];
% a
% % a = [2 4];
% % b = [3 6];
% % d = norm(b-a)
% f = [2 3;
%      4 6;
%      5 12];
%  d = [1 2];
%  f-d
d = [2 4 5 2 3 4;
    67 1 3 1 6 67;
    34 22 17 3 22 55];

dd = sortrows(d',1)