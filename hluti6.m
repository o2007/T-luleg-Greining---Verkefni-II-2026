clear; clc;
disp('');
disp('Hluti 6:');
disp('');

P0 = [0,  0];
P1 = [1,  3];
P2 = [3, -2];
P3 = [4,  1];

[x, y, hradi] = bezierFerill(P0, P1, P2, P3);

quadTol   = 1e-8;
newtonTol = 1e-8;

L = adaptQuad(hradi, 0, 1, quadTol);
fprintf('Heildar bogalengd: %.6f\n\n', L);

teiknaJafnskiptingu(x, y, hradi, L, 6, ...
    'Bezier-ferill: 6 jafnlangir butar', quadTol, newtonTol);

teiknaJafnskiptingu(x, y, hradi, L, 20, ...
    'Bezier-ferill: 20 jafnlangir butar', quadTol, newtonTol);

% Hreyfimynd
fjoldiRamma = 200;
tPlot = linspace(0, 1, 1000);
xFerill = x(tPlot);
yFerill = y(tPlot);

tGildi = linspace(0, 1, fjoldiRamma);
titlar1 = arrayfun(@(t) sprintf('Upprunalegar breytur: t = %.3f', t), ...
                   tGildi, 'UniformOutput', false);
hreyfaFeril(xFerill, yFerill, x(tGildi), y(tGildi), 'r', titlar1);

pause(1);

[xJafnt, yJafnt] = finnaStadsetningar(x, y, hradi, L, fjoldiRamma, quadTol);
sGildi = linspace(0, 1, fjoldiRamma);
titlar2 = arrayfun(@(s) sprintf('Jafnar bogalengdir: s = %.3f', s), ...
                   sGildi, 'UniformOutput', false);
hreyfaFeril(xFerill, yFerill, xJafnt, yJafnt, 'g', titlar2);


function [x, y, hradi] = bezierFerill(P0, P1, P2, P3)
    x = @(t) (1-t).^3 .* P0(1) + 3*(1-t).^2 .* t .* P1(1) + ...
             3*(1-t) .* t.^2 .* P2(1) + t.^3 .* P3(1);

    y = @(t) (1-t).^3 .* P0(2) + 3*(1-t).^2 .* t .* P1(2) + ...
             3*(1-t) .* t.^2 .* P2(2) + t.^3 .* P3(2);

    dx = @(t) 3*(1-t).^2 .* (P1(1)-P0(1)) + ...
              6*(1-t) .* t .* (P2(1)-P1(1)) + ...
              3*t.^2 .* (P3(1)-P2(1));

    dy = @(t) 3*(1-t).^2 .* (P1(2)-P0(2)) + ...
              6*(1-t) .* t .* (P2(2)-P1(2)) + ...
              3*t.^2 .* (P3(2)-P2(2));

    hradi = @(t) sqrt(dx(t).^2 + dy(t).^2);
end


function teiknaJafnskiptingu(x, y, hradi, L, n, titill, quadTol, newtonTol)
    tGildi = jafnskipting(hradi, L, n, quadTol, newtonTol);
    xPunktar = x(tGildi);
    yPunktar = y(tGildi);

    bogaLengdir = zeros(1, n);
    for k = 1:n
        bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
    end

    fprintf('n = %d, hamarks fravik = %.2e\n', ...
            n, max(bogaLengdir) - min(bogaLengdir));

    tPlot = linspace(0, 1, 1000);
    figure;
    plot(x(tPlot), y(tPlot), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(xPunktar, yPunktar, 'bo', 'MarkerSize', 7, ...
         'MarkerFaceColor', 'm');

    for k = 1:n
        plot(xPunktar(k:k+1), yPunktar(k:k+1), 'b-', 'LineWidth', 1);
    end

    for k = 1:n+1
        text(xPunktar(k), yPunktar(k), sprintf('  %d', k-1), 'FontSize', 9);
    end

    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');
    title(titill);
    hold off;
end


function tGildi = jafnskipting(hradi, L, n, quadTol, newtonTol)
    sGildi = linspace(0, 1, n+1);

    tGildi = arrayfun(@(s) Newton( ...
        @(t) adaptQuad(hradi, 0, t, quadTol) - s*L, ...
        @(t) hradi(t), ...
        s, ...
        newtonTol), sGildi);
end


function hreyfaFeril(xFerill, yFerill, xPunktur, yPunktur, litur, titill)
    fjoldiRamma = length(xPunktur);

    figure;
    plot(xFerill, yFerill, 'b-', 'LineWidth', 1.5);
    hold on;
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');

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

    heilt = (b - a) * (f(a) + f(b)) / 2;
    vinstri = (c - a) * (f(a) + f(c)) / 2;
    haegri = (b - c) * (f(c) + f(b)) / 2;

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
        root = root - f(root) / df(root);
    end
end
