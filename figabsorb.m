clear;
mm=200;
nn=200;
fid1=fopen('abr.rec','r');
MA=fread(fid1,[mm,nn],'single');
fclose(fid1);
figure;imagesc(MA);colorbar;

nt=1000;
fid2=fopen('ricker.o','r');
rick=fread(fid2,nt,'single');
fclose(fid2);
figure;plot(rick);