function [ superposition ] = SUP( x ) % �ж��Ƿ������ӵ���̬
    [a, b] = size(x);
    s = x(1).^2;
    for i = 1:a - 1
        s = s + x(i+1).^2;
    end
    if x(1) == 1 && x(a) == 0
        disp('��״̬������̬');
    elseif x(1) == 0 && x(a) == 1
        disp('��״̬�����Ӽ���̬');
    elseif (1 - s) <= (10^-30) && (1-s) >= (-10^-30)
        disp('��״̬�����ӵ���̬');
    else
        disp('��״̬�������ӵ���̬');
    end
end

% SUP(QSS('0'))
% SUP(QSS('1'))
% SUP(0.70/2*QSS('0') + 0.70*QSS('1'))
% SUP(sqrt(2)/2*QSS('0') + sqrt(2)/2*QSS('1'))
% SUP(0.5*QSS('00') + 0.5*QSS('01') + 0.5*QSS('10') + 0.5*QSS('11'))

