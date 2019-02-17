function y = func(x, DIM)
% ------------------------------------------
% Griewank function% The global minima: f(x) = 0
% Position(0)
% ------------------------------------------
n = DIM;
fr = 4000;
s = 0;
p = 1;
for j = 1:n; s = s + x(j)^2; end
for j = 1:n; p = p*cos(x(j)/sqrt(j));end
y = s/fr-p+1;
end

