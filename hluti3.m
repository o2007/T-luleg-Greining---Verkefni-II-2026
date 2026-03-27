%Hluti 3
disp(' ');
disp('Liður 3:');
x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + (0.3 + 1.8*t - 8.1*t.^2).^2);

quadTol = 1e-8;
bisectTol = 1e-8;

L = adaptQuad(hradi, 0, 1, quadTol);
%% n = 4
timiBisect_n4 = keyraFyrirN(4, hradi, L, quadTol, bisectTol, x, y);

%% n = 20
timiBisect_n20 = keyraFyrirN(20, hradi, L, quadTol, bisectTol, x, y);


function timiBisect = keyraFyrirN(n, hradi, L, quadTol, bisectTol, x, y)

    tPlot = linspace(0, 1, 1000);

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

    fprintf('Keyrslutími fyrir n = %d: %.6f s\n', n, timiBisect);

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
    title(sprintf('Ferillinn skiptur í %d jafnlanga búta með Bisection aðferðinni', n));
    hold off;
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
 
