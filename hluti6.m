%% Liður 6
 
P0 = [0,  0];
P1 = [1,  3];
P2 = [3, -2];
P3 = [4,  1];
 
[x, y, hradi] = bezierCurve(P0, P1, P2, P3);
 
quadTol   = 1e-8;
newtonTol = 1e-8;
 
L = adaptQuad(hradi, 0, 1, quadTol);
fprintf('Heildar bogalengd: %.6f\n\n', L);
 
plotEquipartition(x, y, hradi, L, 6,  'Bézier-ferill: 6 jafnlangir bútar',  quadTol, newtonTol);
plotEquipartition(x, y, hradi, L, 20, 'Bézier-ferill: 20 jafnlangir bútar', quadTol, newtonTol);
 
nFrames = 200;
tPlot = linspace(0, 1, 1000);
xFull = x(tPlot);
yFull = y(tPlot);
 
tVals = linspace(0, 1, nFrames);
titlar1 = arrayfun(@(t) sprintf('Upprunulegar breytur: t = %.3f', t), tVals, 'UniformOutput', false);
animateFerill(xFull, yFull, x(tVals), y(tVals), 'r', titlar1);
 
pause(1);
% ATHUGA - Þessi kóði tekur langan tíma að keyra 
[xS, yS] = progressPositions(x, y, hradi, L, @(s) s, nFrames, quadTol);
sVals = linspace(0, 1, nFrames);
titlar2 = arrayfun(@(s) sprintf('Jöfn hraði: s = %.3f', s), sVals, 'UniformOutput', false);
animateFerill(xFull, yFull, xS, yS, 'g', titlar2);

function [x, y, hradi] = bezierCurve(P0, P1, P2, P3)
    x = @(t) (1-t).^3*P0(1) + 3*(1-t).^2.*t*P1(1) + 3*(1-t).*t.^2*P2(1) + t.^3*P3(1);
    y = @(t) (1-t).^3*P0(2) + 3*(1-t).^2.*t*P1(2) + 3*(1-t).*t.^2*P2(2) + t.^3*P3(2);

    dx = @(t) 3*(1-t).^2.*(P1(1)-P0(1)) + 6*(1-t).*t.*(P2(1)-P1(1)) + 3*t.^2.*(P3(1)-P2(1));
    dy = @(t) 3*(1-t).^2.*(P1(2)-P0(2)) + 6*(1-t).*t.*(P2(2)-P1(2)) + 3*t.^2.*(P3(2)-P2(2));

    hradi = @(t) sqrt(dx(t).^2 + dy(t).^2);
end

function plotEquipartition(x, y, hradi, L, n, titill, quadTol, newtonTol)
    tGildi = equipartition(hradi, L, n, quadTol, newtonTol);
    xP = x(tGildi);
    yP = y(tGildi);

    bogaLengdir = zeros(1, n);
    for k = 1:n
        bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
    end
    fprintf('n = %d, hámarks frávik: %.2e\n', n, max(bogaLengdir)-min(bogaLengdir));

    tPlot = linspace(0, 1, 1000);
    figure;
    plot(x(tPlot), y(tPlot), 'r-', 'LineWidth', 1.5); hold on;
    plot(xP, yP, 'bo', 'MarkerSize', 7, 'MarkerFaceColor', 'm');
    for k = 1:n
        plot(xP(k:k+1), yP(k:k+1), 'b-', 'LineWidth', 1);
    end
    for k = 1:n+1
        text(xP(k), yP(k), sprintf('  %d', k-1), 'FontSize', 9);
    end
    axis equal; grid on;
    xlabel('x'); ylabel('y');
    title(titill);
    hold off;
end

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
function I = adaptQuad(f,a,b,tol)
    c = (a+b)/2;

    allt = (b-a)*(f(a)+f(b))/2;
    vinstri  = (c-a)*(f(a)+f(c))/2;
    haegri = (b-c)*(f(c)+f(b))/2;

    if abs(allt - (vinstri + haegri)) < 3*tol
        I = vinstri + haegri;
    else
        I = adaptQuad(f,a,c,tol/2) + adaptQuad(f,c,b,tol/2);
    end
end
function tGildi = equipartition(hradi, L, n, quadTol, newtonTol)
    sGildi = linspace(0, 1, n+1);
    tGildi = arrayfun(@(s) Newton(@(t) adaptQuad(hradi, 0, t, quadTol) - s*L, ...
                                  @(t) hradi(t), s, newtonTol), sGildi);
end
function root = Newton(f, df, x0, TOL)
    root = x0;
    for i = 1:100
        if abs(f(root)) <= TOL
            return;
        end
        root = root - f(root)/df(root);
    end
end

