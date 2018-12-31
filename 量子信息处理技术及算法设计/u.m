function [ unitarity ] = u( x ) % �ж��Ƿ��������
    [m, n] = size(x);
    if m ~= n
        disp('�ξ��󲻾�������');
    else
        I = x*x';
        J = eye(m, n);
        if det(J-I) < (10^-30) && det(J-I) > (-10^-30)
            disp('�ξ����������');
        else
            disp('�ξ��󲻾�������');
        end
    end
end

% I = [1 1; 1 -1]/sqrt(2);
% u(I)