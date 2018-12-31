function [ bloch ] = bl( x ) % Bloch«Ú√Ê
    a = 2 * acos(x(1));
    if a == 0
        b = 0;
    else
        b = acos(real(x(2))/sin(a/2));
    end
    t = linspace(0, pi * 2, 37);
    p = linspace(0, pi, 19);
    [t, p] = meshgrid(t, p);
    r = 1;
    x = r.* sin(p).* cos(t);
    y = r.* sin(p).* sin(t);
    z = r.* cos(p);
    zz = z; zzz = z;
    % zz(p < pi/6) = nan; zzz(p > pi/6) = nan;
    surf(x, y, zz, 'facecolor', [0.8 0.8 0.8], 'edgecolor', 'none');
    hold on
    % surf(x, y, zzz, 'facecolor', [1 1 0], 'edgecolor', 'none');
    camlight; lighting gouraud; alpha .5; axis equal;
    pp = linspace(0, 2 * pi, 19);
    h = linspace(0, 0, 19);
    plot3(sin(pp), cos(pp), h, '--')
    quiver3([0, 0], [0, 0], [0 0 ], [0 sin(a)*cos(b)], [0 sin(a)*sin(b)], [0 cos(a)])
end

% bl(QSS('0'))
% bl(QSS(1))
% c = sqrt(2)/2;
% bl(c*QSS(0) + c*QSS(1));
% bl(c*QSS(0) + exp(i*pi/6)*c*QSS(1));
