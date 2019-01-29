clc
close all

pre_frac = zeros(9, 3, 3);
post_frac = zeros(54, 3, 3);

for j = 1:3
    % Load data from .csv
    pre = clean_csv_read(['output/FFL', num2str(j),'_NGS_pre_counts.csv']);
    post = clean_csv_read(['output/FFL', num2str(j),'_NGS_post_counts.csv']);

    % slice pre data to relevant wells
    pre = pre(1:6:end, :);

    % calculate proportions
    for i = 1:length(pre(:, 1))
        pre_frac(i, :, j) = pre(i, :) ./ sum(pre(i, :));
    end
    for i = 1:length(post(:, 1))
        post_frac(i, :, j) = post(i, :) ./ sum(post(i, :));
    end
end

% define actual ratios
x = [0.1:0.1:0.9]'*(2/3);
y = flipud(x);
z = ones(9,1)*(1/3);
actual_ratios = flipud([x,y,z]);

% plot

for j = 1:3
    figure
for i = 1:9
    subplot(3,3,i)
    line = i;
    location = line*6;
    start = location - 5;
    where = start:location;
    x1 = 1:2;
    x2 = 3:4;
    x3 = 5:6;
    yx = [pre_frac(line, 1, j), mean(post_frac(where, 1, j))];
    yy = [pre_frac(line, 2, j), mean(post_frac(where, 2, j))];
    yz = [pre_frac(line, 3, j), mean(post_frac(where, 3, j))];
    yxa = [actual_ratios(line, 1), actual_ratios(line, 1)];
    yya = [actual_ratios(line, 2), actual_ratios(line, 2)];
    yza = [actual_ratios(line, 3), actual_ratios(line, 3)];
    
    errorx = 2:2:6;
    errory = mean(post_frac(where, :, j));
    errore = std(post_frac(where, :, j));
    
    hold on
    
    plot(x1,yxa,x2,yya,x3,yza,'LineWidth',2,'LineStyle','--','Color','k')
    plot(x1,yx,x2,yy,x3,yz,'LineWidth',2,'Color','k')
    errorbar(errorx, errory, errore, 'LineWidth', 2, 'LineStyle', 'none','Color','k')
    ax = gca;
    ylabel('Fraction')
    tick_locations = [1.5, 3.5, 5.5];
    ax.XTick = tick_locations;
    ax.XTickLabel = {'X', 'Y', 'Z'};
    ylim([0, 0.8])
    box on
    hold off
end
end