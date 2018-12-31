function [tensors] = ts(x, n)  % n次的张量积
    I = x;
	if n > 1 && (n == fix(n))
	    temp = I;
        for i = 1 : n - 1
            temp = t(temp, I);
        end
        I = temp;
        tensors = I;
    else
        disp('请输入整数次')
    end
end

