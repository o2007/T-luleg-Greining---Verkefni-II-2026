clear; clc;

x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;

% Liður 1
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + (0.3 + 1.8*t - 8.1*t.^2).^2);

T = 1;
quadTol = 1e-8;

L = adaptQuad(hradi, 0, T, quadTol);

disp('');
disp('Liður 1:');
disp('');
fprintf('Bogalengd = %.10f\n', L);

% Liður 2
s = 0.5;
bisectTol = 1e-8;

totalLen = adaptQuad(hradi, 0, 1, quadTol);

klukka = tic;
tStjarna = bisectT(hradi, s, totalLen, quadTol, bisectTol);
timiBisect2 = toc(klukka);

disp('');
disp('Liður 2:');
disp('');
fprintf('t* = %.3f\n', tStjarna);

% Liður 3
disp('');
disp('Liður 3:');

% Fyrir n = 4
n = 4;

klukka = tic;
sGildi = linspace(0, 1, n+1);
tGildi = arrayfun(@(s) bisectT(hradi, s, L, quadTol, bisectTol), sGildi);
timiBisect3_n4 = toc(klukka);

xPunktar = x(tGildi);
yPunktar = y(tGildi);

bogaLengdir = zeros(1, n);
for k = 1:n
    bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
end

t3 = linspace(0, 1, 1000);
figure;
plot(x(t3), y(t3), 'r-', 'LineWidth', 1.5);
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
title('Ferillinn skiptur í 4 jafnlanga búta með Bisection aðferðinni');
hold off;

% Fyrir n = 20
n = 20;

klukka = tic;
sGildi = linspace(0, 1, n+1);
tGildi = arrayfun(@(s) bisectT(hradi, s, L, quadTol, bisectTol), sGildi);
timiBisect3 = toc(klukka);

xPunktar = x(tGildi);
yPunktar = y(tGildi);

bogaLengdir = zeros(1, n);
for k = 1:n
    bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
end

t3 = linspace(0, 1, 1000);
figure;
plot(x(t3), y(t3), 'r-', 'LineWidth', 1.5);
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
title('Ferillinn skiptur í 20 jafnlanga búta með Bisection aðferðinni');
hold off;

% Liður 4
disp('');
disp('Liður 4:');
disp('');

newtonTol = 1e-8;
s = 0.5;
t0 = s;

f = @(t) adaptQuad(hradi, 0, t, quadTol) - s*L;
df = @(t) hradi(t);

klukka = tic;
tStjarna = Newton(f, df, t0, newtonTol);
timiNewton2 = toc(klukka);

fprintf('t* = %.3f\n', tStjarna);

% Liður 3 með Newton
n = 4;
klukka = tic;
sGildi = linspace(0, 1, n+1);
tGildi = arrayfun(@(s) Newton(@(t) adaptQuad(hradi, 0, t, quadTol) - s*L, ...
                              @(t) hradi(t), ...
                              s, ...
                              newtonTol), sGildi);
timiNewton3_n4 = toc(klukka);

xPunktar = x(tGildi);
yPunktar = y(tGildi);

bogaLengdir = zeros(1, n);
for k = 1:n
    bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
end

t4 = linspace(0, 1, 1000);
figure;
plot(x(t4), y(t4), 'b-', 'LineWidth', 1.5);
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
title('Ferlinum skipt í 4 jafnlangan búta með Newton aðferðinni');
hold off;

n = 20;

klukka = tic;
sGildi = linspace(0, 1, n+1);
tGildi = arrayfun(@(s) Newton(@(t) adaptQuad(hradi, 0, t, quadTol) - s*L, ...
                              @(t) hradi(t), ...
                              s, ...
                              newtonTol), sGildi);
timiNewton3 = toc(klukka);

xPunktar = x(tGildi);
yPunktar = y(tGildi);

bogaLengdir = zeros(1, n);
for k = 1:n
    bogaLengdir(k) = adaptQuad(hradi, tGildi(k), tGildi(k+1), quadTol);
end

disp('');

t4 = linspace(0, 1, 1000);
figure;
plot(x(t4), y(t4), 'b-', 'LineWidth', 1.5);
hold on;
plot(xPunktar, yPunktar, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
for k = 1:n
    plot(xPunktar(k:k+1), yPunktar(k:k+1), 'r-');
end
axis equal;
grid on;
xlabel('x(t)');
ylabel('y(t)');
title('Ferlinum skipt í 20 jafnlangan búta með Newton aðferðinni');
hold off;

disp('');
disp('Samanburður á keyrslutíma:');
fprintf('Liður 2 Bisection = %.6f s og Newton = %.6f s\n', timiBisect2, timiNewton2);
fprintf('Liður 3 fyrir n=4: Bisection = %.6f s, Newton = %.6f s\n', timiBisect3_n4, timiNewton3_n4);
fprintf('Liður 3 fyrir n=20: Bisection = %.6f s, Newton = %.6f s\n', timiBisect3, timiNewton3);

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
    fa = adaptQuad(f, 0, a, quadTol) - s*totalLen;

    while (b-a)/2 > tol
        c = (a+b)/2;
        fc = adaptQuad(f, 0, c, quadTol) - s*totalLen;

        if fa*fc < 0
            b = c;
        else
            a = c;
            fa = fc;
        end
    end

    tStjarna = (a+b)/2;
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