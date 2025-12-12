clc
clear
[~,basic,element,atom_num,position]=read_POSCAR_car('30_36_zig_relax.vasp');
position=position(:,1)*basic(1,:)+position(:,2)*basic(2,:)+position(:,3)*basic(3,:);
max(position(:,3))
% [counts,centers]=hist(position(:,2))
[dis_y,order_y]=sort(position(:,2));
one_position=position(order_y(1:atom_num/2,:),:);

[dis,order]=sort(one_position(:,1));

new_position=one_position(order,:);
L=sqrt(diff(new_position(:,1)).^2+diff(new_position(:,3)).^2);
x=linspace(0,max(position(:,1)),500);
y=interp1(new_position(:,1),new_position(:,3),x,'spline','extrap');

diff_y=gradient(y)./gradient(x);
diff_yy=gradient(diff_y)./gradient(x);
R=((1+(diff_y).^2).^(3/2))./diff_yy;
k=1./R;
K=movmean(k,15);
max(abs(k))


figure(1)
plot(new_position(:,1),new_position(:,3))
hold on
plot(x,y)
hold off

figure(2)
plot(x,diff_y)
hold on
plot(x,diff_yy)
hold off

figure(3)
plot(x,k,'blue')
hold on
plot(x,K)
hold off