display(1);
aveDiff = zeros(1,3);

cells = load('CB15N_007_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = TwotoCenter(cells);
display(max_min(:,1));
aveDiff(1) = mean(max_min(2,:)-max_min(1,:));



cells = load('BB130_LPho_002.nd2 - s=1 - c=3 - z=0 - t=0_16-Sep-2016_CONTOURS_pill_MESH.mat');
max_min = TwotoCenter(cells);
display(max_min(:,1));
aveDiff(2) = mean(max_min(2,:)-max_min(1,:));


cells = load('LS2821_002_31-Jan-2017_CONTOURS_pill_MESH.mat');
max_min = TwotoCenter(cells);
display(max_min(:,1));
aveDiff(3) = mean(max_min(2,:)-max_min(1,:));



% plot the data
figure
bar_plot = bar(aveDiff);
set(bar_plot(1), 'facecolor', 'k');
% legend('Smaller Pole','Larger Pole');
xlabel('Cell Type'), ylabel('Average Difference Between the Poles');
title('2D Comparison of Pole Pointiness By Two-to-Center');
% end