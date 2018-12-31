function [Hadamards] = HS(x, n)  % n次的Hadamard门
    [a, b] = size(x);
	I = [1 1; 1 -1]/sqrt(2);
	if n > 1 && (n == fix(n))
	    temp = I;
	    for i = 1: n-1
		    temp = t(temp, I);
		end
		I = temp;
		[c, d] = size(I);
		if d == a
		    temp = I * x;
			Hadamards = temp;
	    else
		    error('Hadamard门不能作用此量子态');
		end
	else
	    error('请正确输入参数,作用次数必须大于等于2');
	end
end