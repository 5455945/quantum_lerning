function [tensors] = ts(x, n)  % n�ε�������
    I = x;
	if n > 1 && (n == fix(n))
	    temp = I;
        for i = 1 : n - 1
            temp = t(temp, I);
        end
        I = temp;
        tensors = I;
    else
        disp('������������')
    end
end

