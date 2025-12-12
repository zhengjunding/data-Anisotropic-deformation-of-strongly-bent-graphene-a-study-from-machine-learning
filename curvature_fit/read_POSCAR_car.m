function [lattice,basic,element,atom_num,position]=read_POSCAR_car(name)
fid=fopen(name,'r');
system_name=fgetl(fid);
lattice=fgetl(fid);
basic(1,:)=str2num(fgetl(fid));
basic(2,:)=str2num(fgetl(fid));
basic(3,:)=str2num(fgetl(fid));
element=fgetl(fid);
atom_num=str2num(fgetl(fid));
tmp=fgetl(fid);
% tmp=fgetl(fid);
for i=1:sum(atom_num)
    tmp=fgetl(fid);
    tmp=sscanf(tmp,'%f %f %f');
    position(i,:)=tmp';
end
fclose all;
end
