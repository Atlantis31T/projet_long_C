%% Création d'un jeu de données 

Points=[];

% paramétrage : nombre de points de chaque carré
liste=[0:1/50:1];

for i=1:size(liste,2)
    for j=1:size(liste,2)
        Points=[Points;liste(i),liste(j)];
    end
end


for i=1:size(liste,2)
    for j=1:size(liste,2)
        Points=[Points;liste(i)-2,liste(j)+2];
    end
end


for i=1:size(liste,2)
    for j=1:size(liste,2)
        Points=[Points;liste(i)-4,liste(j)+4];
    end
end


for i=1:size(liste,2)
    for j=1:size(liste,2)
        Points=[Points;liste(i),liste(j)+4];
    end
end



for i=1:size(liste,2)
    for j=1:size(liste,2)
        Points=[Points;liste(i)-4,liste(j)];
    end
end


figure
plot(Points(:,1),Points(:,2),'r.')

save Croix Points
