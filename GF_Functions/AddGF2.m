function sum_out = AddGF2(a1, a2, GF)
%   -a1 and a2 are power form GF(2^m) numbers that will be added together via
%    xor
%   -GF is the enumeration of the specifc field
%   -sum_out is the power form of the determined sum

m = size(GF{1},2);
n = 2^m;
if(a1 < -1 || a1 > n-2)
    error("a1 = %d is not a valid power of alpha in GF(2^%d)\n", a1, m);
    return
end
if(a2 < -1 || a2 > n-2)
    error("a2 = %d is not a valid power of alpha in GF(2^%d)\n", a2, m);
    return
end

add1 = GF{a1+2};
add2 = GF{a2+2};

sum_d = double(xor(add1, add2)); %does sum computation

%finds the equivalent in the GF enumeration
sum_out = cellfun(@(x)isequal(x, sum_d),GF, 'un', 0);
sum_out = find([sum_out{:}] == 1) - 2; %subtracts two to get exp of alpha


end

