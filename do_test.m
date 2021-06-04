function do_test(test_name, R_err, GF, G, C, M, sys_gen)
%DO_TEST Completes a signular decoding test, prints out relevant
%information. May throw an error if decoding algorithm errors. Also
%determines the original message bits that were encoded
%INPUTS:
% test_name - a string containing the name of the test
% R_err     - The received codeword polynomial to decode
% t         - the number of errors the code should correct
% GF        - the array of cells enumerating GF(2^m) symbols as elements
%             in GF(2)[a].
% G         - The generating polynomial used for getting message bits
% C         - The codeword polynomial used for comparing to the corrected
%             codeword (C_hat(x)).
% M         - The original message polynomial
% sys_gen   - Boolean value representing the encoding method (generator or
%             systematic). true = generator method, false = systematic method


fprintf("\n  --------------%s--------------\n",test_name);

fprintf("   The codeword is received corresponding to R(x) shown below.\n");
fprintf("   The codeword R(X) is decoded using the Berlekamp-Massey, Chien, Forney process\n");
fprintf("   The orginal codeword is also shown below for comparison:\n");
if(sys_gen)
    fprintf("   Polynomial forms:\n\t");
    print_poly("C(x)", C, false);
    fprintf("\t");
    print_poly("R(x)", R_err, false);
    fprintf("   Power forms:\n\t");
    print_poly("C[x]", C, true);
    fprintf("\t");
    print_poly("R(x)", R_err, true);
else
    fprintf("   Polynomial forms:\n\t");
    print_poly("C_sys(x)", C, false);
    fprintf("\t");
    print_poly("R(x)", R_err, false);
    fprintf("   Power forms:\n\t");
    print_poly("C_sys[x]", C, true);
    fprintf("\t");
    print_poly("    R(x)", R_err, true);
end

k = size(M,2);
n = size(C,2);
t = (n-k)/2;

try
    [C_hat, failed] = RS_Decoder(R_err, t, GF);
catch decode_error
    failed = true;
end

if(failed)
    error("%s Failed! The following message is produced:\n%s", test_name, decode_error.message);
else
    fprintf("   The received word was decoded!\n");
    fprintf("   Polynomial form:\n\t");
    if(sys_gen)
        print_poly("C_hat(x)", C_hat, false);
        fprintf("   Power form:\n\t");
        print_poly("C_hat[x]", C_hat, true);
        fprintf("   The decoded message is: \n\t");
        %below shouldn't fail because it was already vetted by BerlekampMasseyRS
        [msg_bin, msg] = get_message(C_hat, GF, G);
    else
        print_poly("C_sys_hat(x)", C_hat, false);
        fprintf("   Power form:\n\t");
        print_poly("C_sys_hat[x]", C_hat, true);
        fprintf("   The decoded message is: \n\t");
        %below shouldn't fail because it was already vetted by BerlekampMasseyRS
        [msg_bin, msg] = get_message(C_hat, GF, G, k);
    end
    
    print_poly("Message Bits",msg_bin, true);
    fprintf("\t");
    print_poly("M[x]", msg, true);
    fprintf("  Original Message:\n   ");
    print_poly("Orig_M[x]", M, true);
end

end

