 function [nbcluster] = VisualisationImg(filename, type, filtre)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualisation des images

close all;
% Chargement des données

if (type == 0)
% cas d'une image donnée par un fichier au format image (jpg, ...)
  I1=imread([filename, '.png']);
  I1 = I1(1:50, 1:50);
  I=double(I1);
else
% cas où les données sont déjà stockées au format mat
%load fichier.mat
  name = strcat(filename, '.mat')
  load Taeuber-Arp;
end

disp('nbre de points')
[n1,m1,p] = size(I)

% Nombre de classe à charger

nbclusterSC=importdata('sc5x4/nbclusters');
nbclusterMS=importdata('ms5x4/nbclusters');

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
    IDXSC=zeros(n1*m1,3);
    IDXMS=zeros(n1*m1,3);

    Imag=zeros(n1*m1,3);
    G1=I(:,:,1)';
    G2=I(:,:,2)';
    G3=I(:,:,3)';
    Imag(:,1)=G1(:);
    Imag(:,2)=G2(:);
    Imag(:,3)=G3(:);
    nbgrosclusters = 0;
    for i=1:nbclusterSC
        A = importdata(strcat('sc5x4/cluster.final.',num2str(i)));
        if length(A)>filtre
          nbgrosclusters = nbgrosclusters + 1;
          LevelColor=mean(Imag(A(2:end,1),:));
          IDXSC(A(2:end,1),:)=repmat(floor(LevelColor),size(A,1)-1,1);
        end
    end
    nbgrosclusters = 0;
    for i=1:nbclusterMS
        A = importdata(strcat('ms5x4/cluster.final.',num2str(i)));
        if length(A)>filtre
          nbgrosclusters = nbgrosclusters + 1;
          LevelColor=mean(Imag(A(2:end,1),:));
          IDXMS(A(2:end,1),:)=repmat(floor(LevelColor),size(A,1)-1,1);
        end
    end
    % Visualisation de la comparaison entre image brute et résultat de
    % classification
    figure
    I=uint8(I);
    imshow(I)
    %title(['Original data ' ]);
    saveas(gcf, 'img-original', 'epsc')
    figure
     I4=reshape(IDXSC,m1,n1,3);
     I4=ipermute(I4,[2,1,3]);
     I5=uint8(I4);
     imshow(I5)
     rotate3d on
     %title(['Spectral Clust. 5x4 (' int2str(nbgrosclusters) ')']);
     saveas(gcf, 'img-sc', 'epsc')
     figure
     
     I4=reshape(IDXMS,m1,n1,3);
     I4=ipermute(I4,[2,1,3]);
     I5=uint8(I4);
     imshow(I5)
     rotate3d on
     %title(['Mean-shift 5x4 (' int2str(nbgrosclusters) ')']);
     saveas(gcf, 'img-ms', 'epsc')

    
end



return

