%I=imread('rice.jpg');
load ImageCouleur

dataPts=importdata('ImageCouleur.txt');
[clustCent,data2cluster,cluster2dataCell] = MeanShiftClusterImage(dataPts',100,80);

A=reshape(dataPts(:,3),275,194);
subplot(1,2,1)
imagesc(I1)
T=reshape(data2cluster,275,194);
subplot(1,2,2)
imagesc(T')
%[clustCent,data2cluster,cluster2dataCell] = TestGeom('jain.txt', 8);