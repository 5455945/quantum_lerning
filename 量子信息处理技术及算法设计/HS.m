function [Hadamards] = HS(x, n)  % n�ε�Hadamard��
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
		    error('Hadamard�Ų������ô�����̬');
		end
	else
	    error('����ȷ�������,���ô���������ڵ���2');
	end
end