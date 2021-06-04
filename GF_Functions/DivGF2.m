function quo = DivGF2(dvd, dvs, GF)
%DIVGF2:
%   -dvd is the divided, dvs is the divisor, the input parameters that are to be divided, they are
%   in power form and consist of a single element of GF(2^m)
%   -GF is the enumeration of the GF
%   -quo is the quotient output

m = size(GF{1},2);
n = 2^m;
if(dvd < -1 || dvd > n-2)
    error("dividened = %d is not a valid power of alpha in GF(2^%d)\n", dvd, m);
    return
end
if(dvs < -1 || dvs > n-2)
    error("divisor = %d is not a valid power of alpha in GF(2^%d)\n", dvs, m);
    return
end

if(dvd == -1 || dvs == -1)
    quo = -1;
else
    quo = mod((dvd - dvs) + (n-1),n-1);
end

end

