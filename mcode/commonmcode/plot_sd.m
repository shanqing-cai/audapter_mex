function plot_sd_2(lv, sd, colr)
    for n = 1 : length(lv) - 1
        rectangle('Position', [lv(n), sd(1, n), lv(n + 1) - lv(n), sd(2, n) - sd(1, n)], 'FaceColor', colr, 'EdgeColor', 'none');
        hold on;
    end
return