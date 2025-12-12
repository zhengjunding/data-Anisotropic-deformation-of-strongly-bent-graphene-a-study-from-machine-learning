clc
clear
[lattice,basic,element,atom_num,position]=read_POSCAR_car('graphene_arm.vasp');

%%%%%%%%%%%%%%%%%%%%bending 
num=30;
x=0:0.001:num;
A=1.2;
B=max(x)/2;
C=3;
z=A*exp(-(x-B).^2/2/C^2);
% z=A*sin(2*pi*x./(num)); 
figure(1)
plot(x,z)

%%%%%%%%%%%%%%%%%%%%set rotation center
a=basic(1);b=basic(5);c=basic(9);
center=[a*0.5 b*0.5 0];
%%%%%%%%%%%%%%%%%%%division
rx=a*x; 
rz=c*z-min(c*z);
r=[rx;rz];
j=2;
m=1;
% for i=1:size(x,2)
while  m<size(x,2)
for i=m
        dx=rx-rx(i);
    dz=rz-rz(i);
    dl=sqrt(dx.^2+dz.^2);
    L=abs(dl-a);
    L(1:m-1)=L(1:m-1)+5;
    local=find(L==min(L));
    newr(:,1)=[0;0];
    if min(L)<0.1
        newr(:,j)=r(:,local);
        j=j+1;
   end   
end
     m=local;
     % ttt(j)=local;
end
n=size(newr,2);

for i=1:n-2
   test(i)=sqrt((newr(2,i)-newr(2,i+1))^2+(newr(1,i)-newr(1,i+1))^2);
end


%%%%%%%%%%%%%%%%%%rotation
dx=gradient(newr(1,:));
dz=gradient(newr(2,:));
theta=round(-atand(dz./dx),5);
% n=size(newr,2);
for j=1:sum(atom_num)
for i=1:n-1
rot=[cosd(theta(1,i)) 0 sind(theta(1,i));0 1 0 ;-sind(theta(1,i)) 0 cosd(theta(1,i))];
dat((n-1)*(j-1)+i,:)=rot*([a*position(j,1) b*position(j,2) 0]'-center')+center';
end
end
%%%%%%%%%%%%%%%%%% shift
dat(:,1)=dat(:,1)+repmat(newr(1,1:end-1),1,sum(atom_num))';
dat(:,3)=dat(:,3)+(repmat(newr(2,1:end-1),1,sum(atom_num))'+2);



l=[max(dat(:,1))-min(dat(:,1))+1.4246 0 0;0 b 0;0 0 100];
N=(n-1)*atom_num;

%%%%%%%%%%%%%%%%%%%%%%%%print POSCAR
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
fprintf(fid2,'%f %f %f \n',dat');
sta=fclose(fid2);