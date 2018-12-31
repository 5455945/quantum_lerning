function [ unitarity ] = u( x ) % 判定是否具有酉性
    [m, n] = size(x);
    if m ~= n
        disp('次矩阵不具有酉性');
    else
        I = x*x';
        J = eye(m, n);
        if det(J-I) < (10^-30) && det(J-I) > (-10^-30)
            disp('次矩阵具有酉性');
        else
            disp('次矩阵不具有酉性');
        end
    end
end

% I = [1 1; 1 -1]/sqrt(2);
% u(I)