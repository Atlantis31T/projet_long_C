function [clustCent,data2cluster,cluster2dataCell] = TestGeom(filename, bandwidth)

A=importdata(filename);
data=A(2:end,:);
[clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(data',bandwidth);

numClust = length(cluster2dataCell);
x = data';

figure(10),clf,hold on
cVec = 'bgrcmyk';%, cVec = [cVec cVec];
nbC = length(cVec);
nbgros = 0;
for k = 1:numClust
    j = mod(k, nbC);
    myMembers = cluster2dataCell{k};
    if (length(myMembers) > 20)
        nbgros = nbgros + 1;
        myClustCen = clustCent(:,k)
        plot(x(1,myMembers),x(2,myMembers),[cVec(j+1) '.']);
        plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(j+1), 'MarkerSize',10);
    end
end
title(['no shifting, numClust:' int2str(nbgros)])

