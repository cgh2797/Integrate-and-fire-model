function prettyfigure()
% enhances the appearance of the figure plots

ax = gca;
fg = gcf;
set(ax, 'Box', 'Off');
set(ax, 'LineWidth',2);
set(ax, 'FontSize', 14);
set(fg, 'Color', 'w');