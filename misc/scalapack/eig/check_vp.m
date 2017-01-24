clear all;
format long;

Asp = load('Asp3.txt');
vp = load('vp3.txt');

n = 356;
A = reshape(Asp, n, n);

norm(A- A')


vpar = reshape(vp(:,2), n, n);
vseq = reshape(vp(:,1), n, n);

[V,D] = eig(A);

for i = 1:n

    j = n-i+1;

    err_mat(i) = norm(A*V(:,j) - D(j,j)*V(:,j));

    err_par(i) = norm(A*vpar(:,i) - D(j,j)*vpar(:,i));
    err_seq(i) = norm(A*vseq(:,i) - D(j,j)*vseq(:,i));
end

norm(err_mat)
norm(err_seq)
norm(err_par)