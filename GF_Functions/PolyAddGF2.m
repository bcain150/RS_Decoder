function sum = PolyAddGF2(p1,p2, GF)
%POLYADDGF2 adds two polynomial functions in power form outputs the sum
% of those polynomial functions. Add occurs element by element wise
% matching with respective coef locations
% INPUTS:
%   p1, p2 - two polynomial functions in polar form
%   GF - The cell array mapping of groups of GF(2) elements to GF(2^m)
%   elements
% OUPUTS:
%   sum - the power form representation of the polynomial sum


p1_size = size(p1, 2);
p2_size = size(p2, 2);

%note: L is the amount of elements that should be summed together. The
%rest of the elements of the larger polynomial are kept in the sum

if (p1_size > p2_size)
    sum = p1;
    L = p2_size;
else
    sum = p2;
    L = p1_size;
end

for i = 0:L-1
    sum(1,end-i) = AddGF2(p1(1,end-i), p2(1,end-i), GF);
end

end

