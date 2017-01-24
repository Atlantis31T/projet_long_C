 function [nbcluster] = VisualisationImg(basename, type, filtre, method)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des images
%% basename : basename du fichier contenant les données
%% type : type des données originales : 0 == png, 1 == mat
%% filtre : taille minimum des clusters
%% method : 1 == Spectral Clustering, 2 == Mean shift 

close all;
% Chargement des données

if (type == 0)
% cas d'une image donnée par un fichier au format image (jpg, ...)
  I1=imread([basename, '.png']);
  %I1 = I1(1:50, 1:50);
  I=double(I1);
else
% cas où les données sont déjà stockées au format mat
%load fichier.mat -> variable I
  load(basename);
end

disp('nbre de points')
[n1,m1,p] = size(I)

% Nombre de classe à charger
disp('nombre de classes')
nbcluster=importdata('nbclusters');


if p==1
    disp('image niveau de gris')
   % Récupération des clusters et définition du niveau de gris moyen par
   % classe
    IDX=zeros(n1*m1,1); 
    I1=I';
    Imag=I1(:);

    for i=1:nbcluster
        A=importdata(strcat('cluster.final.',num2str(i)));
        LevelColor=mean(Imag(A(2:end,1),:));
        IDX(A(2:end,1))=repmat(LevelColor,size(A,1)-1,1);
    end

  
     % Visualisation de la comparaison entre image brute et résultat de
     % classification
    figure
    hold on
    subplot(1,2,1)
    I=uint8(I);
    imshow(I)
    title(['Original data']);
    
    subplot(1,2,2)
    I4=reshape(IDX,m1,n1);
    I4=ipermute(I4,[2,1,3]);
    I5=uint8(I4);
    imshow(I5)
    rotate3d on
    title(['Clusters number = ' int2str(nbcluster)]);
     
   
else
    
    disp('image couleur')
   % Récupération des clusters et définition de la couleur par
   % classe
    IDX=zeros(n1*m1,3);
    Imag=zeros(n1*m1,3);
    G1=I(:,:,1)';
    G2=I(:,:,2)';
    G3=I(:,:,3)';
    Imag(:,1)=G1(:);
    Imag(:,2)=G2(:);
    Imag(:,3)=G3(:);
        nbgrosclusters = 0;
    for i=1:nbcluster
        A=importdata(strcat('cluster.final.',num2str(i)));
        if length(A)>filtre
          nbgrosclusters = nbgrosclusters + 1;
          LevelColor=mean(Imag(A(2:end,1),:));
          IDX(A(2:end,1),:)=repmat(floor(LevelColor),size(A,1)-1,1);
        end
    end
    
    % Visualisation de la comparaison entre image brute et résultat de
    % classification
    figure
    hold on
    subplot(1,2,1)
    I=uint8(I);
    imshow(I)
     title(['Original data ' ]);
%     saveas(gcf, 'img-original', 'epsc')
     subplot(1,2,2)
     I4=reshape(IDX,m1,n1,3);
     I4=ipermute(I4,[2,1,3]);
     I5=uint8(I4);
     imshow(I5)
     rotate3d on
    if method == 1
      title(['Spectral Clustering 5x4 (' int2str(nbgrosclusters) ')']);
      saveas(gcf, 'img-sc', 'epsc')
    else
      title(['Mean-shift 5x4 (' int2str(nbgrosclusters) ')']);
      saveas(gcf, 'img-ms', 'epsc')
    end
    
end



return

