nbproc = [ 18; 32; 45; 60; 64 ];

nbpoints = [22000;12500;9000;6800;6300];

par_time = [36616;7243;2808;1153;1030];

tot_time = [38927;9415;5127;3495;3157];

speedup = par_time(1) ./ par_time;
speedupt = tot_time(1) ./ tot_time;

eff = (speedup .* nbproc)/18;
cout = par_time ./ nbpoints;
coutt = tot_time ./nbpoints;

ccout  = par_time ./ (881*441);

figure(1)
plot(nbproc, cout, 'r.-','LineWidth', 2, 'MarkerSize',15);
hold on;
plot(nbproc, coutt, 'k-.','LineWidth', 2, 'MarkerSize',15);
legend('parallel cost', 'total cost');
xlabel('number of processors', 'FontSize', 14);
ylabel('computational cost of one point', 'FontSize', 14);
saveas(1, 'cost.eps', 'psc2');

figure(2)
plot(nbproc, speedupt, 'r-','LineWidth', 2, 'MarkerSize',15);
xlabel('number of processors', 'FontSize', 14);
ylabel('speedup (reference time: 18 processors)', 'FontSize', 14);
saveas(2, 'speedup.eps', 'psc2');

figure(3)
plot(nbproc, eff, 'k.-','LineWidth', 2, 'MarkerSize',15);
xlabel('number of processors', 'FontSize', 14);
ylabel('efficiency (reference time: 18 processors)', 'FontSize', 14);
saveas(2, 'efficiency.eps', 'eps');

figure(4)
plot(nbproc, ccout, 'k.-','LineWidth', 2, 'MarkerSize',15);
xlabel('number of processors', 'FontSize', 14);
ylabel('cout', 'FontSize', 14);
saveas(2, 'cooo.eps', 'eps');

