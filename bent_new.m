clc
clear
[lattice,basic,element,atom_num,position]=read_POSCAR_car('graphene.vasp');

%%%%%%%%%%%%%%%%%%%%beding曲线
num=11;
x=0:0.001:num;
A=0.3;
B=0;
C=3;
z=A*exp(-(x-B).^2/2/C^2); %%高斯
% z=A*sin(2*pi*x./(num));  %%三角函数
figure(1)
plot(x,z)
%%%%%%%%%%%%%%%%%%%%线段长度接近整数
% for i=1:size(x,2)-1
%    q(i)=sqrt((x(i)-x(i+1))^2+(z(i)-z(i+1))^2);
% end
% disp(sum(q))
%%%%%%%%%%%%%%%%%%%%
a=basic(1);b=basic(5);c=basic(9);
center=[a*0.5 b*0.5 0];
%%%%%%%%%%%%%%%%%%%寻找临近点
rx=a*x;
rz=c*z-min(c*z);
r=[rx;rz];
% j=2;
% m=1;
% % for i=1:size(x,2)
% while  m<size(x,2)
% for i=m
%     dx=rx-rx(i);
%     dz=rz-rz(i);
%     dl=sqrt(dx.^2+dz.^2);
%     local=knnsearch(dl',a);
%     % L=abs(dl-a);
%     % L(1:m-1)=L(1:m-1)+5;
%     % local=find(L==min(L));
%     % newr(:,1)=[0;0];
%     % if min(L)<0.1
%         
%         j=j+1;
%    % end   
% end
%      m=local;
% end

% dl=sqrt(rx.^2+rz.^2);
tmp=[diff(rx) ;diff(rz)];
dl=[0 cumsum(sqrt(sum(tmp.^2)))];
N=max(dl)/a;
disp(N)
local=a.*(1:fix(N)-1);
local_position=knnsearch(dl',local');
newr=r(:,local_position);

% dl=zeros(1,length(rx));
% for i=1:length(rx)-1
%     dl(i)=
% end
%


n=size(newr,2);
% figure(2)
% scatter(newr(1,:),newr(2,:))
% hold on
% plot(a*x,c*z)
% hold off
% % equal axis
% for i=1:n-2
%    test(i)=sqrt((newr(2,i)-newr(2,i+1))^2+(newr(1,i)-newr(1,i+1))^2);
% end
% figure(5)
% plot(test)

%%%%%%%%%%%%%%%%%%转动
dx=gradient(newr(1,:));
dz=gradient(newr(2,:));
theta=round(-atand(dz./dx),5);
% n=size(newr,2);
for j=1:sum(atom_num)
for i=1:n
rot=[cosd(theta(1,i)) 0 sind(theta(1,i));0 1 0 ;-sind(theta(1,i)) 0 cosd(theta(1,i))];
dat((n)*(j-1)+i,:)=rot*([a*position(j,1) b*position(j,2) 0]'-center')+center';
end
end

dat(:,1)=dat(:,1)+repmat(newr(1,1:end),1,sum(atom_num))'-a;
dat(:,3)=dat(:,3)+(repmat(newr(2,1:end),1,sum(atom_num))'+2);
dat(:,1)=dat(:,1)-min(dat(:,1));

inver_dat=[-dat(:,1) dat(:,2) dat(:,3) ];
inver_dat(find(inver_dat(:,1)==0),:)=[];

all_dat=[dat ;inver_dat];
all_dat(:,1)=all_dat(:,1)-min(all_dat(:,1));

l=[2*(num-1)*a-2*a 0 0;0 b 0;0 0 100];
N=length(all_dat);
disp(num-1)
disp(N/4)

%%%%%%%%%%%%%%%%%%%%%%%%写出POSCAR
outfile=strcat('POSCAR','.vasp');
fid2=fopen(outfile,'w');
% firstLine
fprintf(fid2,'bending_stru');
fprintf(fid2,'\n');
fprintf(fid2,'%f \n',1);
fprintf(fid2,'%f %f %f \n',l);
fprintf(fid2,element);
fprintf(fid2,'\n');
fprintf(fid2,'%d \n',N);
fprintf(fid2,'C \n');
fprintf(fid2,'%f %f %f \n',all_dat');
sta=fclose(fid2);