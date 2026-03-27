clear; clc;

% Hraði ferilsins
hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + ...
                  (0.3 + 1.8*t - 8.1*t.^2).^2);

quadTol = 1e-8;
T = 1;

L = adaptQuad(hradi, 0, T, quadTol);

disp('');
disp('Hluti 1:');
disp('');
fprintf('Bogalengd = %.10f\n', L);


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
