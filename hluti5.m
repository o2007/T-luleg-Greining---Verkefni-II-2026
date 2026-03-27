% Hluti 5 -- Hreyfimynd
clear; clc;
disp('');
disp('Liður 5:');

x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + ...
                  (0.3 + 1.8*t - 8.1*t.^2).^2);

quadTol = 1e-8;
L = adaptQuad(hradi, 0, 1, quadTol);

fjoldiRamma = 200;
tPlot = linspace(0, 1, 1000);
xFerill = x(tPlot);
yFerill = y(tPlot);

% Hreyfing með jofnum skrefum i t
tGildi = linspace(0, 1, fjoldiRamma);
titlar1 = arrayfun(@(t) sprintf('Upprunalegar breytur: t = %.3f', t), ...
                   tGildi, 'UniformOutput', false);

hreyfaFeril(xFerill, yFerill, x(tGildi), y(tGildi), 'r', titlar1);

pause(1);

% Hreyfing med jofnum skrefum i bogalengd
[xJafnt, yJafnt] = finnaStadsetningar(x, y, hradi, L, fjoldiRamma, quadTol);

sGildi = linspace(0, 1, fjoldiRamma);
titlar2 = arrayfun(@(s) sprintf('Jafnur hraði: s = %.3f', s), ...
                   sGildi, 'UniformOutput', false);

hreyfaFeril(xFerill, yFerill, xJafnt, yJafnt, 'g', titlar2);


function hreyfaFeril(xFerill, yFerill, xPunktur, yPunktur, litur, titill)
    fjoldiRamma = length(xPunktur);

    figure;
    plot(xFerill, yFerill, 'b-', 'LineWidth', 1.5);
    hold on;
    axis equal;
    grid on;
    xlabel('x(t)');
    ylabel('y(t)');

    punktur = plot(NaN, NaN, 'o', 'Color', litur, ...
                   'MarkerFaceColor', litur, 'MarkerSize', 8);

    for k = 1:fjoldiRamma
        set(punktur, 'XData', xPunktur(k), 'YData', yPunktur(k));
        if iscell(titill)
            title(titill{k});
        else
            title(titill);
        end
        drawnow;
        pause(0.01);
    end
end


function [xPunktur, yPunktur, tStjarna] = finnaStadsetningar(x, y, hradi, L, fjoldiRamma, quadTol)
    sGildi = linspace(0, 1, fjoldiRamma);
    tStjarna = zeros(1, fjoldiRamma);

    for k = 1:fjoldiRamma
        if sGildi(k) <= 0
            tStjarna(k) = 0;
        elseif sGildi(k) >= 1
            tStjarna(k) = 1;
        else
            tStjarna(k) = Newton(@(t) adaptQuad(hradi, 0, t, quadTol) - sGildi(k)*L, ...
                                 @(t) hradi(t), sGildi(k), 1e-8);
        end
    end

    xPunktur = x(tStjarna);
    yPunktur = y(tStjarna);
end


function I = adaptQuad(f, a, b, tol)
    c = (a + b)/2;

    heilt = (b - a)*(f(a) + f(b))/2;
    vinstri = (c - a)*(f(a) + f(c))/2;
    haegri = (b - c)*(f(c) + f(b))/2;

    if abs(heilt - (vinstri + haegri)) < 3*tol
        I = vinstri + haegri;
    else
        I = adaptQuad(f, a, c, tol/2) + adaptQuad(f, c, b, tol/2);
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
