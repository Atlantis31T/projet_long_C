 function [nbcluster] = VisualisationGeom(filtre, sizepoint)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des figures géométriques
close all

% Chargement des données
Data0=importdata(['r15/R15.txt']);
DataR15=Data0(2:end, :);
Data0=importdata(['d31/D31.txt']);
DataD31=Data0(2:end, :);
% Taille des données
disp('nbre de points')
p = 2;

nbclusters(1) = importdata('r15/sc2x2/nbclusters');
nbclusters(2) = importdata('r15/ms2x2/nbclusters');
nbclusters(3) = importdata('d31/sc2x2/nbclusters');
nbclusters(4) = importdata('d31/ms2x2/nbclusters');
nbclusters
nbcluster=max(nbclusters)

if p==2
   disp('Figure 2D')
   cm=colormap(hsv(10*(nbcluster+1)));
   
   % Visualisation résultats sc 1x1
   nbcluster = importdata('r15/sc2x2/nbclusters');
   subplot(1,4,1)
   
   %axis ([0 70 25 45]);
   axis square
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('r15/sc2x2/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(DataR15(ii(2:end),1),DataR15(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
   title(['Spectral Clust. 2x2 (' int2str(nbgrosclusters) ')']);
   ylabel('R15', 'FontSize', 12, 'FontWeight', 'bold');
   
  % Visualisation résultats sc 2x2
   nbcluster = importdata('r15/ms2x2/nbclusters');
   subplot(1,4,2) 
   %axis ([0 70 25 45]);
   axis square
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('r15/ms2x2/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(DataR15(ii(2:end),1),DataR15(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
   title(['Mean-Shift 2x2 (' int2str(nbgrosclusters) ')']);
   
   % Visualisation résultats ms 1x1
   nbcluster = importdata('d31/sc2x2/nbclusters');
   subplot(1,4,3) 
   %axis ([0 70 25 45]);
   axis square
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('d31/sc2x2/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(DataD31(ii(2:end),1),DataD31(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
  
   %set(h,'position',[x y])
   h = title(['Spectral Clust. 2x2 (' int2str(nbgrosclusters) ')']);
   ylabel('D31', 'FontSize', 12, 'FontWeight', 'bold');
   get(h,'position')
    
  % Visualisation résultats ms 2x2
   nbcluster = importdata('d31/ms2x2/nbclusters');
   subplot(1,4,4) 
   %axis ([0 70 25 45]);
   axis square
   
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('d31/ms2x2/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(DataD31(ii(2:end),1),DataD31(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
   title(['Mean-Shift 2x2 (' int2str(nbgrosclusters) ')']);
   
   saveas(gcf, 'R15D31', 'epsc')
    
else
    disp('Figure 3D')
    cm=colormap(hsv(10*(nbcluster+1)));
    % Visualisation des données brutes  
    subplot(1,2,1)    
    hold on 
    plot3(Data(:,1),Data(:,2),Data(:,3),'g.','MarkerSize',15);
     title(['Original data ' ]);
   % Récupération et Visualisation des clusters 
    subplot(1,2,2)
        hold on 
    for i=1:nbcluster
        ii=importdata(strcat('cluster.final.',num2str(i)));
        plot3(Data(ii(2:end),1),Data(ii(2:end),2),Data(ii(2:end),3),'.','color',[cm(i*10,:)],'MarkerSize',15);
    end
     title(['Clusters number = ' int2str(nbcluster)]);
    
end

return
