function [controlnot] = CNOT(x)   % �������ܿط���
    y = x;
	[a, b] = size(x);
	if b == 2
	    if (y(1) == '0' || y(1) == '1') && (y(2) == '0' || y(2) == '1')  % �����ж�
		    c = num2str(xor(str2num(y(1)), str2num(y(2))));
			controlnot = [y(1) c];
		else
		    error('��������λ0��1�ַ����');
		end
    else
	    error('�ܿط��Ų������ô�����̬');
	end
end