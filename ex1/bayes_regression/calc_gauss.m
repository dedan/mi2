function [P] = calc_gauss(w,m,sig)

if size(w,2) < size(m,2)
   w = repmat(w,[1 size(m,2)]);
elseif size(m,2) < size(w,2)
   m = repmat(m,[1 size(w,2)]);
end

N = size(w,1);

P = [];
for i = 1:size(w,2)
  P(i) = 1/((2*pi)^(N/2) * sqrt(det(sig))) * exp( -1/2*(w(:,i)-m(:,i))'*inv(sig)*(w(:,i)-m(:,i)));
end

if N > 1
 P = reshape(P,repmat(length(P)^(1/N),[1 N]));
end;
