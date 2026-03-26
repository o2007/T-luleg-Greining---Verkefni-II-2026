clear; clc; close all;

hradi = @(t) sqrt((0.3 + 7.8*t - 14.1*t.^2).^2 + (0.3 + 1.8*t - 8.1*t.^2).^2);

quadTol = 1e-8;
bisectTol = 1e-8;
s = 0.5;

totalLen = adaptQuad(hradi, 0, 1, quadTol);

klukka = tic;
tStjarna = bisectT(hradi, s, totalLen, quadTol, bisectTol);
timiBisect = toc(klukka);

disp(' ');
disp('Liður 2:');
disp(' ');
fprintf('t* = %.3f\n', tStjarna);
fprintf('Keyrslutími = %.6f s\n', timiBisect);


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