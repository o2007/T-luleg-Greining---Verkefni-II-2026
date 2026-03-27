clear; clc;

% Skilgreina ferilinn
x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;

% Hraði ferilsins
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + ...
                  (0.3 + 1.8*t - 8.1*t.^2).^2);

% Vikmörk
T = 1;
quadTol = 1e-8;
bisectTol = 1e-8;
newtonTol = 1e-8;

%% Hluti 1
L = adaptQuad(hradi, 0, T, quadTol);

disp('');
disp('Hluti 1:');
disp('');
fprintf('Bogalengd = %.10f\n', L);

%% Hluti 2
s = 0.5;

klukka = tic;
tStjarna = bisectT(hradi, s, L, quadTol, bisectTol);
timiBisect2 = toc(klukka);

disp('');
disp('Hluti 2:');
disp('');
fprintf('t* = %.3f\n', tStjarna);

%% Hluti 3
disp('');
disp('Hluti 3:');
disp('');

timiBisect3_n4 = hluti3Bisect(4, hradi, L, quadTol, bisectTol, x, y);
timiBisect3_n20 = hluti3Bisect(20, hradi, L, quadTol, bisectTol, x, y);

%% Hluti 4
disp('');
disp('Hluti 4:');
disp('');

f = @(t) adaptQuad(hradi, 0, t, quadTol) - s*L;
df = @(t) hradi(t);

klukka = tic;
tStjarnaNewton = Newton(f, df, s, newtonTol);
timiNewton2 = toc(klukka);

fprintf('t* = %.3f\n', tStjarnaNewton);

timiNewton_n4 = keyraNewton(4, hradi, L, quadTol, newtonTol, x, y);
timiNewton_n20 = keyraNewton(20, hradi, L, quadTol, newtonTol, x, y);

disp('');
disp('Samanburdur a keyrslutima:');
fprintf('Hluti 2 Bisection = %.6f s og Newton = %.6f s\n', timiBisect2, timiNewton2);
fprintf('Hluti 3 fyrir n = 4: Bisection = %.6f s, Newton = %.6f s\n', timiBisect3_n4, timiNewton_n4);
fprintf('Hluti 3 fyrir n = 20: Bisection = %.6f s, Newton = %.6f s\n', timiBisect3_n20, timiNewton_n20);


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


function tStjarna = bisectT(f, s, totalLen, quadTol, tol)
    if s == 0
        tStjarna = 0;
        return
    elseif s == 1
        tStjarna = 1;
        return
    end

    a = 0;
    b = 1;
    fa = -s * totalLen;

    while (b - a)/2 > tol
        c = (a + b)/2;
        fc = adaptQuad(f, 0, c, quadTol) - s * totalLen;

        if fa * fc < 0
            b = c;
        else
            a = c;
            fa = fc;
        end
    end

    tStjarna = (a + b)/2;
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


function timiBisect = hluti3Bisect(n, hradi, L, quadTol, bisectTol, x, y)

    klukka = tic;
    sGildi = linspace(0, 1, n+1);
    tGildi = arrayfun(@(s) bisectT(hradi, s, L, quadTol, bisectTol), sGildi);
    timiBisect = toc(klukka);

    xPunktar = x(tGildi);
    yPunktar = y(tGildi);

    bogaLengdir = zeros(1, n);
    for k = 1:n
        bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
    end

    fprintf('Keyrslutimi fyrir n = %d: %.6f s\n', n, timiBisect);

    tPlot = linspace(0, 1, 1000);
    figure;
    plot(x(tPlot), y(tPlot), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(xPunktar, yPunktar, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'm');

    for k = 1:n
        plot(xPunktar(k:k+1), yPunktar(k:k+1), 'b-');
    end

    for k = 1:n+1
        text(xPunktar(k), yPunktar(k), sprintf('  %d', k-1), 'FontSize', 9);
    end

    axis equal;
    grid on;
    xlabel('x(t)');
    ylabel('y(t)');
    title(sprintf('Ferillinn skiptur i %d jafnlanga buta med Bisection adferdinni', n));
    hold off;
end


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

    bogaLengdir = zeros(1, n);
    for k = 1:n
        bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
    end

    fprintf('Keyrslutimi fyrir n = %d med Newton: %.6f s\n', n, timiNewton);

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
    title(sprintf('Ferillinn skiptur i %d jafnlanga buta med Newton adferdinni', n));
    hold off;
end
