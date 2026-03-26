% Liður 5 — Hreyfimynd
clear; clc;
disp('');
disp('Liður 5:');
 
x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + ...
                  (0.3 + 1.8*t - 8.1*t.^2).^2);
T = 1;
quadTol = 1e-8;
L = adaptQuad(hradi, 0, T, quadTol);
nFrames = 200;
tPlot = linspace(0, 1, 1000);
xFull = x(tPlot);
yFull = y(tPlot);
 
tVals = linspace(0, 1, nFrames);
titlar1 = arrayfun(@(t) sprintf('Upprunalegar breytur: t = %.3f', t), tVals, 'UniformOutput', false);
animateFerill(xFull, yFull, x(tVals), y(tVals), 'r', titlar1);
 
pause(1);
 
[xS, yS] = progressPositions(x, y, hradi, L, @(s) s, nFrames, quadTol);
sVals = linspace(0, 1, nFrames);
titlar2 = arrayfun(@(s) sprintf('Jöfn hraði: s = %.3f', s), sVals, 'UniformOutput', false);
animateFerill(xFull, yFull, xS, yS, 'g', titlar2);

function animateFerill(xFull, yFull, xBall, yBall, litur, titill)
    nFrames = length(xBall);
    xlims = [min(xFull)-0.5, max(xFull)+0.5];
    ylims = [min(yFull)-0.5, max(yFull)+0.5];

    figure;
    plot(xFull, yFull, 'b-', 'LineWidth', 1.5);
    set(gca, 'XLim', xlims, 'YLim', ylims, 'Drawmode', 'fast', 'Visible', 'on');
    axis square;

    ball = line('color', litur, 'Marker', 'o', 'MarkerSize', 10, ...
                'LineWidth', 2, 'erase', 'xor', 'xdata', [], 'ydata', []);

    for k = 1:nFrames
        set(ball, 'xdata', xBall(k), 'ydata', yBall(k));
        if iscell(titill)
            title(titill{k});
        else
            title(titill);
        end
        drawnow;
        pause(0.01);
    end
end


function [xBall, yBall, tStars] = progressPositions(x, y, hradi, L, C, nFrames, quadTol)
    sVals = linspace(0, 1, nFrames);
    cVals = C(sVals);
    tStars = zeros(1, nFrames);

    for k = 1:nFrames
        if cVals(k) <= 0
            tStars(k) = 0;
        elseif cVals(k) >= 1
            tStars(k) = 1;
        else
            tStars(k) = Newton(@(t) adaptQuad(hradi, 0, t, quadTol) - cVals(k)*L, ...
                               @(t) hradi(t), cVals(k), 1e-8);
        end
    end

    xBall = x(tStars);
    yBall = y(tStars);
end
