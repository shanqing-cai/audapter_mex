function plot_sd_1(af1, af2, sdf2, colr)
% Assume ascending af2.
    for n = 1 : length(af1) - 1
%         if (sdf1(n+1)+af1(n+1)>af1(n)-sdf1(n))
            rectangle('Position', [af1(n), af2(n)-sdf2(n),af1(n+1)-af1(n),sdf2(n+1)+af2(n+1)-(af2(n)-sdf2(n))],...
                'FaceColor', colr, 'EdgeColor', 'none');
%         else
%             rectangle('Position', [sdf1(n+1)+af1(n+1),af2(n),(af1(n)-sdf1(n))-(sdf1(n+1)+af1(n+1)),af2(n+1)-af2(n)],...
%                 'FaceColor', colr, 'EdgeColor', 'none');
%         end
%         hold on;
    end
return