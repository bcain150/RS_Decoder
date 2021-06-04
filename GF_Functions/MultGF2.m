function prod_out= MultGF2(m1, m2, GF)
%MULTGF2:
%   -m1 and m2 are the input parameters that are to be multiplied, they are
%   in power form and consist of a single element of GF(2^m)
%   -GF is the enumeration of the GF

m = size(GF{1},2);
n = 2^m;
if(m1 < -1 || m1 > n-2)
    error("a1 = %d is not a valid power of alpha in GF(2^%d)\n", a1, m);
    return
end
if(m2 < -1 || m2 > n-2)
    error("a2 = %d is not a valid power of alpha in GF(2^%d)\n", a2, m);
    return
end

if(m1 == -1 || m2 == -1)
    prod_out = -1;
else
    prod_out = mod(m1 + m2, n-1);
end

end

