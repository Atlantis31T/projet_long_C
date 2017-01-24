 function [nbcluster] = VisualisationGeom(filename, filtre)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des figures géométriques
close all

% Chargement des données
Data=importdata([filename '.txt']);

% Taille des données
disp('nbre de points')
[n1,p]=size(Data)

% Nombre de classe à charger
disp('nombre de classes')
nbcluster = importdata('nbclusters');



if p==2
   disp('Figure 2D')
   cm=colormap(hsv(10*(nbcluster+1)));
   % Visualisation des données brutes
    subplot(1,2,1)  
    hold on
    plot(Data(:,1),Data(:,2),'k.','MarkerSize',20);
     title(['Original data ' ]);
   % Récupération et Visualisation des clusters 
    subplot(1,2,2)
    hold on
    nbgrosclusters = 0;
    for i=1:nbcluster
        ii=importdata(strcat('cluster.final.',num2str(i)));
        if(length(ii(2:end)) > filtre)
           nbgrosclusters = nbgrosclusters + 1;
           plot(Data(ii(2:end),1),Data(ii(2:end),2),'.','color',[cm(i*10,:)],'MarkerSize',20);
           pause
        end
        %pause
    end
     title(['Clusters number = ' int2str(nbgrosclusters)]);
    
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
