% TwotoCenter
% cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
% cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
% cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');

max_min = TwotoCenter(cells);

% perform a paired t-test on the radii values at the two poles
[h,p,ci,stats] = ttest(max_min(2,:),max_min(1,:));

display(h);
display(p);
display(ci);
display(stats);

ymin = min(max_min(:)) - 0.5;
ymax = max(max_min(:)) + 0.5;


% plot the data
figure
bar_plot = bar(transpose(max_min));
% xlim([1,100]);
ylim([ymin, ymax]);
set(bar_plot(1), 'facecolor', 'k');
set(bar_plot(2), 'facecolor', 'r');
legend('Smaller Pole','Larger Pole');
xlabel('Cell #'), ylabel('proxy_width');
title('2D Comparison of Pole Pointiness By Two-to-Center');

%-------------------------------------------------------------------------%
% TwoPointWidth

% cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
% cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
% cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');

max_min = TwoPointWidth(cells);
% perform a paired t-test on the radii values at the two poles
[h,p,ci,stats] = ttest(max_min(2,:),max_min(1,:));

display(h);
display(p);
display(ci);
display(stats);

ymin = min(max_min(:)) - 0.5;
ymax = max(max_min(:)) + 0.5;

% for i = 0:int8(numCells/50)-2
% plot the data
figure
bar_plot = bar(transpose(max_min));
% xlim([1,100]);
ylim([ymin, ymax]);
set(bar_plot(1), 'facecolor', 'k');
set(bar_plot(2), 'facecolor', 'r');
legend('Smaller Pole','Larger Pole');
xlabel('Cell #'), ylabel('proxy_width');
title('2D Comparison of Pole Pointiness by 2 point width approximation');

%-------------------------------------------------------------------------%
% AreaPerimRatio (larger value is LESS pointy!)
cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
% cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
% cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');

max_min = AreaPerimRatio(cells);
% perform a paired t-test on the radii values at the two poles
[h,p,ci,stats] = ttest(max_min(2,:),max_min(1,:));

display(h);
display(p);
display(ci);
display(stats);

ymin = 0;
ymax = max(max_min(:)) + 0.2;

% for i = 0:int8(numCells/50)-2
% plot the data
figure
bar_plot = bar(transpose(max_min));
% xlim([1,100]);
ylim([ymin, ymax]);
set(bar_plot(1), 'facecolor', 'k');
set(bar_plot(2), 'facecolor', 'r');
legend('Smaller Pole','Larger Pole');
xlabel('Cell #'), ylabel('area/perim (greater = less pointy)');
title('2D Comparison of Pole Pointiness by area/perim ratio');





