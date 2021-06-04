function [C_hat, failure] = RS_Decoder(R_GF2m, t, GF, prnt_flag)
%{
RS_DECODER: implementation of the Berlekamp-Massey algorithm used
for decoding and correcting a received Reed-Solomon codeword in GF(2^m).
This function calls other functions to derive these polynomials needed for
error detection and correction: the syndrome polynomial = S(x) 
(get_syndromes.m), the error locator polynomial = lambda(x)
(get_error_loc.m), the error location polynomial = the inverse roots of 
lambda(x) (chien_search.m), the error polynomial = e(x) (coeff given by 
evalForney.m). These functions utilize the GF(2^m) arithmetic functions 
created in proj3. Additional GF(2^m) functions were created when a 
specific computation is needed for a given algorithm.

INPUTS:
 R_GF2m    - the received codeword polynomial in power form to be decoded
 t         - an integer value denoting the max amount of errors the code is
             able to correct.
 GF        - GF is the array of cells enumerating GF(2^m) symbols as
             elements in GF(2)[a].
 prnt_flag - optional Boolean parameter that defaults to false, displays 
             algorithm results step-by-step if true, display is omitted if 
             false.
OUTPUTS:
 C_hat   - The corrected codeword in power format that is returned at the 
           end of the function.
 failure - Boolean output. If true, decoder failed and C_hat is set to
           equal all -1's because a^-1 represents a^inf = 0 (C_hat becomes
           the all zeros codeword)
%}

%see if need to set default value of prnt_flag
if ~exist('prnt_flag','var')
    prnt_flag = false;
end

%set initial conditions
failure = false;
n = size(R_GF2m,2);
m = size(GF{1},2);
k = n-2*t;
if(prnt_flag)
    fprintf("---------Reed-Solomon Decoding---------\n");
    fprintf("\nPART 1 - Initial Conditions:\n");
    fprintf("  1.) Decoding a (%d,%d) %d-error correcting Reed-Solomon Code\n", n, k, t);
    fprintf("    with operations over GF(2^%d)\n", m);
    fprintf("  2.) The received codeword in both power and polynomial form:\n");
    print_poly("    R[x]", R_GF2m, true);
    print_poly("    R(x)", R_GF2m, false);
    fprintf("  3.) The array of cells enumerating GF(2^%d) symbols as\n", m);
    fprintf("    elements in GF(2)[a] is shown below:\n");
    print_poly("    a^inf", GF{1}, true);
    for elm=2:2^m
        fprintf("    a^%d", elm-2);
        print_poly("  ", GF{elm}, true);
    end
end

%get syndromes
if(prnt_flag)
    fprintf("\nPART 2 - Construct the syndrome polynomial S(x):\n");
end
s = get_syndromes(t, R_GF2m, GF, prnt_flag);

%if syndromes are all inf then there are no errors, the codeword passes,
%and the algorithm ends here
if(s(:) == -1)
    C_hat = R_GF2m;
    if(prnt_flag)
        fprintf("  2.) After the computation of all syndrome coefficients, S(x) is constructed\n");
        fprintf("    in power form as: ");
        print_poly("S[x]", s, true);
        fprintf("  3.) All S(x) coefficients are a^inf and thus no errors were detected\n");
        fprintf("  4.) Therefore, the corrected codeword C_hat is set to R(x)\n");
        fprintf("    and the algorithm ends:\n");
        print_poly("    C_hat[x]", C_hat, true);
        print_poly("    C_hat(x)", C_hat, false);
    end
    
    return;
end

if(prnt_flag)
    fprintf("  2.) After the computation of all syndrome coefficients, S(x) is constructed:\n");
    print_poly("    S[x]", s, true);
    print_poly("    S(x)", s, false);
    fprintf("  3.) Since some syndrome coefficients are non-zero, there are errors in R(x).\n");
    fprintf("\nPART 3 - Use the Berlekamp-Massey Algrorithm to compute lambda(x):\n");
end

%get the error location polynomial

lambda = get_error_loc(s, GF, prnt_flag);
deg_lambda = size(lambda,2) - 1;

%if the degree of lambda is greater than t then the decoder fails and the
%algorithm ends here.
if(deg_lambda > t)
    failure = true;
    C_hat = zeros(1,n);
    C_hat(:) = -1;
    
    if(prnt_flag)
        fprintf("  3.) The error locator polynomial lambda(x) for R(x):\n");
        print_poly("   lambda[x]", lambda, true);
        print_poly("   lambda(x)", lambda, false);
        fprintf("  4.) However the degree of lambda(x)(%d) is greater than t(%d)\n", deg_lambda, t);
        fprintf("   and so a decoder failure is triggered.");
        fprintf("  5.) C_hat is set to the all zeros code word (a^inf) and the algorithm quits\n");
        print_poly("   C_hat[x]", C_hat, true);
    else
        error("The degree of lambda(x)(%d) is greater than t(%d)! Decoder Failed!\n", deg_lambda, t);
    end
    
    return;
end

if(prnt_flag)
    fprintf("  3.) The error locator polynomial lambda(x) for R(x):\n");
    print_poly("   lambda[x]", lambda, true);
    print_poly("   lambda(x)", lambda, false);
    fprintf("\nPART 5 - Determine the error locations of R(x) using Chien Search and lambda(x):\n");
end

%get the error locations via chien search
try
    %could error
    error_locations = chien_search(lambda, GF, prnt_flag);
catch chien_error
    %if the number of errors is less than the degree of lambda then the 
    %decoder fails and the algorithm ends here.
    failure = true;
    C_hat = zeros(1,n);
    C_hat(:) = -1;
    
    if(prnt_flag)
        fprintf("  2.) The Chien Search algorithm abruptly stopped with the error message:\n");
        fprintf("    %s", chien_error.message);
        fprintf("  3.) C_hat is set to the all zeros code word (a^inf) and the algorithm quits\n");
        print_poly(" C_hat[x]", C_hat, true);
    else
        err_msg = sprintf("The Chien Search algorithm abruptly stopped with the error message:\n%sDecoder Failed!\n",chien_error.message);
        error("%s", err_msg);
    end
    return;
end

num_errors = size(error_locations, 2);


s_plus_1 = s;
s_plus_1(end) = 0;
%error evaluator polynomial below...
omega = PolyMultGF2(lambda, s_plus_1, GF);
omega = omega(1,end-2*t:end); %do this to take mod(x^2t+1)

%need derivative of error locator (lambda_prime)
lambda_prime = FormalDerivGF2(lambda, GF);

if(prnt_flag)
    fprintf("\nPART 6 - Determine the error polynomial e(x):\n");
    fprintf(" - In order to determine the error polynomial, first determine\n");
    fprintf("   the error evaluator polynomial, omega(x) using the Forney method.\n");
    fprintf(" - Also formally derive lambda(x) to get lambda'(x) for its\n");
    fprintf("   use in the Forney equation.\n");
    fprintf("  1.) Omega is given by the equation: omega(x) = lambda(x)*(S(x) + 1)mod(x^(r+1))\n");
    print_poly("    omega[x]", omega, true);
    print_poly("    omega(x)", omega, false);
    fprintf("  2.) lambda'(x) is computed and shown below:\n");
    print_poly("    lambda'[x]", lambda_prime, true);
    print_poly("    lambda'(x)", lambda_prime, false);
    fprintf("  3.) Use the Forney equation of the form:\n");
    fprintf("    e_i = -X_i*omega(X_i^-1)/lambda'(X_i^-1), where\n");
    fprintf("    e_i is the error value for error i, and X_i is a root of \n");
    fprintf("    the error locator polynomial.\n");
    fprintf("  4.) Evaluate the Forney equation for all error locations:\n");
end

e = zeros(1, n); %initialize size of the error polynomial
e(:) = -1; %fill with a^inf

for e_loc = 1:num_errors
    X = error_locations(1, e_loc);
    e_idx = n - X;
    e(1, e_idx) = evalForney(omega, lambda_prime, X, GF);
    if(prnt_flag)
        fprintf("    for error %d at x^%d, the error value is: %d\n", e_loc, X, e(1, e_idx));
    end
end

C_hat = PolyAddGF2(e, R_GF2m, GF);


if(prnt_flag)
    fprintf("  5.) Construct the the error polynomial, e(x), using the error values\n");
    fprintf("    and the previously determined error locations:\n");
    print_poly("    e[x]", e, true);
    print_poly("    e(x)", e, false);
    fprintf("\nPART 7 - Finally, add e(x) to the received codeword polynomial R(x):\n");
    fprintf(" - This gives the corrected codeword polynomial, C_hat(x), shown below:\n");
    print_poly("    C_hat[x]", C_hat, true);
    print_poly("    C_hat(x)", C_hat, false);
end

end

