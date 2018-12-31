function [controlnot] = CNOT(x)   % 多量子受控非门
    y = x;
	[a, b] = size(x);
	if b == 2
	    if (y(1) == '0' || y(1) == '1') && (y(2) == '0' || y(2) == '1')  % 输入判定
		    c = num2str(xor(str2num(y(1)), str2num(y(2))));
			controlnot = [y(1) c];
		else
		    error('请输入两位0、1字符组合');
		end
    else
	    error('受控非门不能作用此量子态');
	end
end