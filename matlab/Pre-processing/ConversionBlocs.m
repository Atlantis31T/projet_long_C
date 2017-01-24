load Croix.mat

fid=fopen('Croix.txt','w')

n=size(Points,1);
m=size(Points,2);

fprintf(fid, '%i %i \n', n,m);

for i=1:n
        fprintf(fid, '%12.8f %12.8f \n' ,Points(i,1),Points(i,2));
end

fclose(fid);
disp('fini');