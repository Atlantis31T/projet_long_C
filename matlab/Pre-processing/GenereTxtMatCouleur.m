function [nbpoints] = GenereTxtMatCouleur(filename)
% Génération des datas pour SC parallèle

I = imread([filename '.jpg']); % exemple d'image à charger
%I = double(I);
fid=fopen([filename '.txt'],'w'); % fichier data à créer;


n = size(I,1);
m = size(I,2);
fprintf(fid, '%i %i \n', 2, 3); % matrice 2D = en couleur
fprintf(fid, '%i %i \n', n,m);  % dimensions de la matrice

for i=1:n
    for j=1:m
        fprintf(fid, '%12.8f %12.8f %12.8f\n', I(i,j,1), I(i,j,2), I(i,j,3));
    end
end
fclose(fid);
disp('fini');

save(filename,'I');

nbpoints = n*m;