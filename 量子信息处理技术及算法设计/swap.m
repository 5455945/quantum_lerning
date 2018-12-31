function [ swap ] = swap( x, m, n)
    [a, b] = size(x);
    if (b < n) || (b < m) || (m == n) || n ~= fix(n) || m ~= fix(m)
        disp('请正确输入要交换的位')
    else
        temp = x(m);
        x(m) = x(n);
        x(n) = temp;
        swap = x;
    end
end

% swap('0111', 1, 2)
% swap('0111', 2, 1)
% swap('0111', 2, 1.2)
% swap('0111', 2, 5)