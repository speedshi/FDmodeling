clear;
m=101;% number of traces
n=1000;% number of samples
fid1=fopen('vz.rec','r');
A1=fread(fid1,[m,n],'single');
fid2=fopen('ompvz.rec','r');
A2=fread(fid2,[m,n],'single');
fclose(fid1);
fclose(fid2);
A1=A1';
figure;wigb(A1);
figure;imagesc(A1);colorbar;
A2=A2';
figure;wigb(A2);
figure;imagesc(A2);colorbar;
figure;wigb(A1-A2);
figure;imagesc(A1-A2);colorbar;
ntrace=50;% trace number which you want to extract
xt=1:1:n;
ya1=A1(:,ntrace); 
ya2=A2(:,ntrace);
figure;plot(ya1);
figure;plot(ya2);
figure;h1=plot(xt,ya1,'r',xt,ya2,'b');
xlabel('Time(ms)');ylabel('VZ amplitude(normalized)');
legend('se','omp');
set(h1,'linewidth',1.8);