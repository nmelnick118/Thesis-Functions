% TwotoCenter

diffs = cell(1,3);
display(diffs);
aveDiff = zeros(1,3);

cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = TwotoCenter(cells);
diffs{1} = max_min(2,:)-max_min(1,:);
aveDiff(1) = mean(diffs{1});

cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
max_min = TwotoCenter(cells);
diffs{2} = max_min(2,:)-max_min(1,:);
aveDiff(2) = mean(diffs{2});

cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = TwotoCenter(cells);
diffs{3} = max_min(2,:)-max_min(1,:);
aveDiff(3) = mean(diffs{3});

ymax = max(aveDiff) + 0.1;

% plot the data
figure
bar_plot = bar(aveDiff);
set(bar_plot(1), 'facecolor', 'k');
ylim([0, ymax]);
% legend('Smaller Pole','Larger Pole');
str = {'CB15N', 'BB130', 'LS2821'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
ylabel('Average Difference Between the Poles');
title('2D Comparison of Pole Pointiness By Two-to-Center');

% statistical significance
[h,p,ci,stats] = ttest2(diffs{1}, diffs{2});
display(h);
%display(p);
%display(ci);
%display(stats);

[h,p,ci,stats] = ttest2(diffs{2}, diffs{3});
display(h);
%display(p);
%display(ci);
%display(stats);

[h,p,ci,stats] = ttest2(diffs{3}, diffs{1});
display(h);
%display(p);
%display(ci);
%display(stats);

%-------------------------------------------------------------------------%
% TwoPointWidth

diffs = cell(1,3);
display(diffs);
aveDiff = zeros(1,3);

cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = TwoPointWidth(cells);
diffs{1} = max_min(2,:)-max_min(1,:);
aveDiff(1) = mean(diffs{1});

cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
max_min = TwoPointWidth(cells);
diffs{2} = max_min(2,:)-max_min(1,:);
aveDiff(2) = mean(diffs{2});

cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = TwoPointWidth(cells);
diffs{3} = max_min(2,:)-max_min(1,:);
aveDiff(3) = mean(diffs{3});

ymax = max(aveDiff) + 0.1;

% plot the data
figure
bar_plot = bar(aveDiff);
set(bar_plot(1), 'facecolor', 'k');
ylim([0, ymax]);
% legend('Smaller Pole','Larger Pole');
str = {'CB15N', 'BB130', 'LS2821'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
ylabel('Average Difference Between the Poles');
title('2D Comparison of Pole Pointiness By Two-Point-Width');


% statistical significance
[h,p,ci,stats] = ttest2(diffs{1}, diffs{2});
display(h);
%display(p);
%display(ci);
%display(stats);

[h,p,ci,stats] = ttest2(diffs{2}, diffs{3});
display(h);
%display(p);
%display(ci);
%display(stats);

[h,p,ci,stats] = ttest2(diffs{3}, diffs{1});
display(h);
%display(p);
%display(ci);
%display(stats);

%-------------------------------------------------------------------------%
% AreaPerimRatio

diffs = cell(1,3);
display(diffs);
aveDiff = zeros(1,3);

cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = AreaPerimRatio(cells);
diffs{1} = max_min(2,:)-max_min(1,:);
aveDiff(1) = mean(diffs{1});

cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
max_min = AreaPerimRatio(cells);
diffs{2} = max_min(2,:)-max_min(1,:);
aveDiff(2) = mean(diffs{2});

cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = AreaPerimRatio(cells);
diffs{3} = max_min(2,:)-max_min(1,:);
aveDiff(3) = mean(diffs{3});

ymax = max(aveDiff) + 0.1;

% plot the data
figure
bar_plot = bar(aveDiff);
set(bar_plot(1), 'facecolor', 'k');
ylim([0, ymax]);
% legend('Smaller Pole','Larger Pole');
str = {'CB15N', 'BB130', 'LS2821'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
ylabel('Average Difference Between the Poles');
title('2D Comparison of Pole Pointiness by area/perim ratio');


% statistical significance

[h,p,ci,stats] = ttest2(diffs{1}, diffs{2});
display(h);
%display(p);
%display(ci);
%display(stats);

[h,p,ci,stats] = ttest2(diffs{2}, diffs{3});
display(h);
%display(p);
%display(ci);
%display(stats);

[h,p,ci,stats] = ttest2(diffs{3}, diffs{1});
display(h);
%display(p);
%display(ci);
%display(stats);




