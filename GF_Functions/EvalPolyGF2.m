function output = EvalPolyGF2(poly, x, GF)
%EVALPOLYGF2
%   -poly is a power form array polynomial and x is the power form input to
%    that polynomial
%   -GF is the enumeration of the given GF
%   -output is the power form output of the function given the input x
%   -Evaluates the polynomial at x using MultGF2 and AddGF2 functions

m = size(GF{1},2);
n = 2^m;
deg = size(poly,2) - 1;
for i=1:deg
    if(x == -1 || poly(1,i) == -1)
        poly(1,i) = -1;
    else
        temp = mod((deg+1-i)*x,(n-1));
        poly(1,i) = MultGF2(temp, poly(1,i),GF);
    end
end
output = poly(1,1);

for i=2:deg+1
    output = AddGF2(output,poly(1,i),GF);
end
    

end

