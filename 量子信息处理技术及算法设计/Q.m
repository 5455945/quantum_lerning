function [ quantum] = Q( x )
    i = x;
    if i == 0
        quantum = [1;0];
    else
        if i == 1
            quantum = [0;1];
        else
            disp('x ����Ϊ0��1');
        end
    end
end

