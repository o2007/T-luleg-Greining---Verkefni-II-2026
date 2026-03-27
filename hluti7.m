x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;
 
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + ...
                  (0.3 + 1.8*t - 8.1*t.^2).^2);
 
quadTol = 1e-8;
nFrames = 200;
 
L = adaptQuad(hradi, 0, 1, quadTol);
 
tPlot = linspace(0, 1, 1000);
xFull = x(tPlot);
yFull = y(tPlot);
sVals = linspace(0, 1, nFrames);
 
% Progress-kúrfur og nöfn
progressKurfur = { @(s) s,                              'C(s) = s (jafn hraði)';
                   @(s) s.^(1/3),                       'C(s) = s^{1/3}';
                   @(s) s.^2,                           'C(s) = s^2';
                   @(s) sin(s*pi/2),                    'C(s) = sin(s\pi/2)';
                   @(s) 0.5 + 0.5*sin((2*s-1)*pi/2),   'C(s) = ease-in-out' };
 
for i = 1:size(progressKurfur, 1)
    C    = progressKurfur{i, 1};
    nafn = progressKurfur{i, 2};
 
    [xBall, yBall] = progressPositions(x, y, hradi, L, C, nFrames, quadTol);
 
    cVals = C(sVals);
    titlar = arrayfun(@(k) sprintf('%s:  s = %.3f,  C(s) = %.3f', nafn, sVals(k), cVals(k)), ...
                      1:nFrames, 'UniformOutput', false);
 
    animateFerill(xFull, yFull, xBall, yBall, 'r', titlar);
    pause(1);
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
function root = Newton(f, df, x0, TOL)
    root = x0;
    for i = 1:100
        if abs(f(root)) <= TOL
            return;
        end
        root = root - f(root)/df(root);
    end
end
