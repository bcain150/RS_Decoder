function derivative = FormalDerivGF2(poly, GF)
%{
FORMALDERIVGF2 Computes the formal derivative of a power form polynomial
in GF(2^m)
INPUTS:
    poly - the power form polynomial in GF(2^m)
    GF - Array of cells mapping the binary groups to GF(2^m) alpha coeff
OUTPUT:
    derivative - the formal derivative of poly
%}

deg_poly = size(poly,2) - 1;
%power form array to store derivative (one element shorter than poly
%because derivative takes away one degree)
derivative = zeros(1,deg_poly); 
x_pwr = deg_poly; %used to keep track of current position x^i element in poly

for i = 1:deg_poly
    %if the current power in poly is even then the coeff cancel in GF(2^m)
    if(mod(x_pwr, 2) == 0)
        derivative(i) = -1;
    else %otherwise it's an odd multiplier to the coeff and all cancel except 1
        derivative(i) = poly(i);
    end
    x_pwr = x_pwr - 1;
end

%Before return happens remove padding -1's or infs
for i=1:deg_poly
    if(derivative(1,i) ~= -1)
        break;
    end
end

derivative = derivative(1,i:end);

end

