 function [nbcluster] = VisualisationGeom(filename, filtre, sizepoint)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des figures géométriques
close all

% Chargement des données
Data0=importdata([filename '.txt']);
Data=Data0(2:end, :);

% Taille des données
disp('nbre de points')
[n1,p]=size(Data)

nbclusters(1) = importdata('sc1/nbclusters');
nbclusters(2) = importdata('sc2x2/nbclusters');
nbclusters(3) = importdata('ms1/nbclusters');
nbclusters(4) = importdata('ms2x2/nbclusters');
nbclusters
nbcluster=max(nbclusters)

if p==2
   disp('Figure 2D')
   cm=colormap(hsv(10*(nbcluster+1)));
   
   % Visualisation résultats sc 1x1
   nbcluster = importdata('sc1/nbclusters')
   subplot(1,4,1)
   
   %axis ([0 70 25 45]);
   axis square
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('sc1/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(Data(ii(2:end),1),Data(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
   title(['Spectral Clust. 1x1 (' int2str(nbgrosclusters) ')']);
   ylabel(filename, 'FontSize', 12, 'FontWeight', 'bold');
   
  % Visualisation résultats sc 2x2
   nbcluster = importdata('sc2x2/nbclusters');
   subplot(1,4,2) 
   %axis ([0 70 25 45]);
   axis square
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('sc2x2/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(Data(ii(2:end),1),Data(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
   title(['Spectral Clust. 2x2 (' int2str(nbgrosclusters) ')']);
   
   % Visualisation résultats ms 1x1
   nbcluster = importdata('ms1/nbclusters');
   subplot(1,4,3) 
   %axis ([0 70 25 45]);
   axis square
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('ms1/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(Data(ii(2:end),1),Data(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
  
   %set(h,'position',[x y])
   h = title(['Mean-Shift 1x1 (' int2str(nbgrosclusters) ')']);
   get(h,'position')
    
  % Visualisation résultats ms 2x2
   nbcluster = importdata('ms2x2/nbclusters');
   subplot(1,4,4) 
   %axis ([0 70 25 45]);
   axis square
   
   hold on
   nbgrosclusters = 0;
   for i=1:nbcluster
       ii=importdata(strcat('ms2x2/cluster.final.',num2str(i)));
       if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;           
           plot(Data(ii(2:end),1),Data(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',sizepoint);
       end
   end
   title(['Mean-Shift 2x2 (' int2str(nbgrosclusters) ')']);
   
   saveas(gcf, filename, 'epsc')
    
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
