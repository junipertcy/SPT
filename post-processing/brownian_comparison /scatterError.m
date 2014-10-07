% scatterError.m (Downloaded from the web)
%
% Function that (more elegantly) draws the scatter plots with
% errorbars at x-y dimension
%
%
%
%
%
% Tzu-Chi Yen
% last modified Sep 2014
%
% originally written by
% Brandon Barker 01/20/2014

function scatterError(x, y, xe, ye, varargin)

nD = length(x);

%Make these defaults later:
dotColor = [1 0.3 0.3]; % conservative pink
yeColor = [0, 0.4, 0.8]; % bright navy blue
xeColor = [0.35, 0.35, 0.35]; % not-too-dark grey

dotSize = 23;

figure();
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gca, 'FontSize', 23);
hold all;

for i = 1:nD
    plot([(x(i) - xe(i)) (x(i) + xe(i))], [y(i) y(i)], 'Color', xeColor);
    %plot([x(i) x(i)], [(y(i) - ye(i)) (y(i) + ye(i))], 'Color', yeColor);
end

scatter(x, y, dotSize, repmat(dotColor, nD, 1));
set(gca, varargin{:});
axis square;