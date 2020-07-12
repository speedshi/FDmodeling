clear;
n1=200;% n1 and n2 depend on the kind of the section. If it's a x-y section, then n1 is the number of grid in Y direction; n2 is the number of grid in X direction;
n2=200;% If it's a y-z section, then n1 is the number of grid in Z direction; n2 is the number of grid in Y direction. Others the same. Note the order is Z-Y-X in the WU model.
n3=100;% n3 is the number of time samples
dname='vzpz.snap';
fid1=fopen(dname,'r');
A1=zeros(n1,n2,n3);
for i=1:n3
     A1(:,:,i)=fread(fid1,[n1,n2],'single');
end
fclose(fid1);

dname2='ompvzpz.snap';
fid2=fopen(dname2,'r');
A2=zeros(n1,n2,n3);
for i=1:n3
     A2(:,:,i)=fread(fid2,[n1,n2],'single');
end
fclose(fid2);

max(max(max(abs(A1-A2))))

% movie(mov,1);
% movie2avi(mov,'vxsnapmx100.avi');
