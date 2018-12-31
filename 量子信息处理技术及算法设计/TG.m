function [ Toffoli ] = TG( x ) % Toffoli门
    y = x;
    [a, b] = size(x);
    if b == 3
        if (y(1) == '0' || y(1) == '1') && (y(2) == '0' || y(2) == '1') && (y(3) == '0' || y(3) == '1')
            c = num2str(xor(str2num(y(3)), and (str2num(y(1)), str2num(y(2)))));
            Toffoli = [y(1) y(2) c];
        else
            disp('请输入三位0、1字符组合')
        end
    else
        disp('受控非门不能作用此量子态')
    end
end

% a = ['000'; '001'; '010'; '011'; '100'; '101'; '110'; '111'];
% for i = 1: 8
%     b = TG(a(i, :));
%     sprintf('%c %c %c 通过TG门后的值为 %c %c %c \n', a(i, 1), a(i, 2), a(i, 3),
%     b(1), b(2), b(3))
% end
