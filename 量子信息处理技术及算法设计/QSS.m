function quantumss = QSS(x) % 多量子比特
y = x;
[a, b] = size(x);
temp = 0;
for i = 1:b
    if (y(i) ~= '0') && (y(i) ~= '1')
        temp = 1;
    end
end

if temp == 0
    I = Q(str2num(y(1)));
    if b == 1
        quantumss = I;
    else
        for i = 1:b-1
            I = t(I, Q(str2num(y(i+1))));
        end
        quantumss = I;
    end
else
    SS = class(x);
    if SS(:, 1:4) =='doub'
        x = num2str(x);
        y = x;
        [a, b] = size(x);
        flag = 0;
        for i = 1:b
            if (y(i) ~= '0') && (y(i) ~= '1')
                flag = 1;
            end
        end
        if flag == 0
            I = Q(str2num(y(1)));
            if b == 1
                quantumss = I;
            else
                for i = 1:b-1
                    I = t(I, Q(str2num(y(i+1))));
                end
                quantumss = I;
            end
        else
            temp = str2num(x);
            temp = dec2bin(temp);
            [c, d] = size(temp);
            I = Q(str2num(temp(1)));
            for i = 1:d-1
                I = t(I, Q(str2num(temp(i+1))));
            end
            quantumss = I;
        end
    else
        temp = str2num(x);
        temp = dec2bin(temp);
        [c, d] = size(temp);
        I = Q(str2num(temp(1)));
        for i = 1:d-1
            I = t(I, Q(str2num(temp(i+1))));
        end
        quantumss = I;
    end
end
end

