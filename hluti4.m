%Hluti 4
clear; clc; close all;

x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + (0.3 + 1.8*t - 8.1*t.^2).^2);

quadTol = 1e-8;
newtonTol = 1e-8;

L = adaptQuad(hradi, 0, 1, quadTol);

disp(' ');
disp('Liður 4:');
disp(' ');

% Liður 2 með Newton
s = 0.5;
t0 = s;

f = @(t) adaptQuad(hradi, 0, t, quadTol) - s*L;
df = @(t) hradi(t);

klukka = tic;
tStjarna = Newton(f, df, t0, newtonTol);
timiNewton2 = toc(klukka);

fprintf('t* = %.3f\n', tStjarna);
fprintf('Keyrslutími fyrir lið 2 með Newton = %.6f s\n', timiNewton2);
%% n = 4
timiNewton_n4 = keyraNewton(4, hradi, L, quadTol, newtonTol, x, y);

%% n = 20
timiNewton_n20 = keyraNewton(20, hradi, L, quadTol, newtonTol, x, y);

fprintf('Liður 2 með Newton = %.6f s\n', timiNewton2);
fprintf('Liður 3 fyrir n = 4 með Newton = %.6f s\n', timiNewton_n4);
fprintf('Liður 3 fyrir n = 20 með Newton = %.6f s\n', timiNewton_n20);


function timiNewton = keyraNewton(n, hradi, L, quadTol, newtonTol, x, y)

    klukka = tic;
    sGildi = linspace(0, 1, n+1);

    tGildi = arrayfun(@(s) Newton( ...
        @(t) adaptQuad(hradi, 0, t, quadTol) - s*L, ...
        @(t) hradi(t), ...
        s, ...
        newtonTol), sGildi);

    timiNewton = toc(klukka);

    xPunktar = x(tGildi);
    yPunktar = y(tGildi);

    tPlot = linspace(0, 1, 1000);

    figure;
    plot(x(tPlot), y(tPlot), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(xPunktar, yPunktar, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

    for k = 1:n
        plot(xPunktar(k:k+1), yPunktar(k:k+1), 'r-');
    end

    for k = 1:n+1
        text(xPunktar(k), yPunktar(k), sprintf('  %d', k-1), 'FontSize', 9);
    end

    axis equal;
    grid on;
    xlabel('x(t)');
    ylabel('y(t)');
    title(sprintf('Ferillinn skiptur í %d jafnlanga búta með Newton aðferðinni', n));
    hold off;

end

function root = Newton(f, df, x0, TOL)
    root = x0;
    for i = 1:100
        if abs(f(root)) <= TOL
            return
        end
        root = root - f(root)/df(root);
    end
end


function I = adaptQuad(f, a, b, tol)
    c = (a+b)/2;

    allt = (b-a)*(f(a)+f(b))/2;
    vinstri = (c-a)*(f(a)+f(c))/2;
    haegri = (b-c)*(f(c)+f(b))/2;

    if abs(allt - (vinstri + haegri)) < 3*tol
        I = vinstri + haegri;
    else
        I = adaptQuad(f, a, c, tol/2) + adaptQuad(f, c, b, tol/2);
    end
end
