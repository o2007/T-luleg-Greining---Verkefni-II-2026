clear; clc;
disp('');
disp('Hluti 7:');
disp('');

x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;

hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + ...
                  (0.3 + 1.8*t - 8.1*t.^2).^2);

quadTol = 1e-8;
fjoldiRamma = 200;

L = adaptQuad(hradi, 0, 1, quadTol);

tPlot = linspace(0, 1, 1000);
xFerill = x(tPlot);
yFerill = y(tPlot);
sGildi = linspace(0, 1, fjoldiRamma);

stodufoll = {
    @(s) s,                            'C(s) = s (jafn hradi)';
    @(s) s.^(1/3),                     'C(s) = s^(1/3)';
    @(s) s.^2,                         'C(s) = s^2';
    @(s) sin(s*pi/2),                  'C(s) = sin(pi*s/2)';
    @(s) 0.5 + 0.5*sin((2*s-1)*pi/2), 'C(s) = ease-in-out'
};

for i = 1:size(stodufoll, 1)
    C = stodufoll{i, 1};
    nafn = stodufoll{i, 2};

    [xPunktur, yPunktur] = finnaStadsetningarC( ...
        x, y, hradi, L, C, fjoldiRamma, quadTol);

    cGildi = C(sGildi);
    titlar = arrayfun(@(k) sprintf('%s: s = %.3f, C(s) = %.3f', ...
                      nafn, sGildi(k), cGildi(k)), ...
                      1:fjoldiRamma, 'UniformOutput', false);

    hreyfaFeril(xFerill, yFerill, xPunktur, yPunktur, 'r', titlar);
    pause(1);
end


function [xPunktur, yPunktur, tStjarna] = finnaStadsetningarC( ...
    x, y, hradi, L, C, fjoldiRamma, quadTol)

    sGildi = linspace(0, 1, fjoldiRamma);
    cGildi = C(sGildi);
    tStjarna = zeros(1, fjoldiRamma);

    for k = 1:fjoldiRamma
        if cGildi(k) <= 0
            tStjarna(k) = 0;
        elseif cGildi(k) >= 1
            tStjarna(k) = 1;
        else
            tStjarna(k) = Newton( ...
                @(t) adaptQuad(hradi, 0, t, quadTol) - cGildi(k)*L, ...
                @(t) hradi(t), cGildi(k), 1e-8);
        end
    end

    xPunktur = x(tStjarna);
    yPunktur = y(tStjarna);
end


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
