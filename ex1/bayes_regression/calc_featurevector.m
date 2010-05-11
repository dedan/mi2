function [phi] = calc_featurevector(f,x)

N = length(f);
for n = 1:N
  phi(n,:) = f{n}(x);
end 
