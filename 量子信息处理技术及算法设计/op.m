function [outerproduct] = op(x, y)
    [m, n] = size(x);
	[p, q] = size(y);
    X = sprintf('%d %d %d %d', m, n, p, q);
    disp(X);
	if n ~= 1 || q ~= 1 || m ~= p
	    error('����ȷ�����������2������')
    else
	    outerproduct = x*y';
    end
end
