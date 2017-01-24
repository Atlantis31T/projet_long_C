N = 4;

A =  hilb(N) + diag([1:-1/N:1/N]);

v = [];
for j = 1:N
v = [v A(j,:)];
end
v=v';
save('A4.txt', 'v', '-ascii', '-double')

