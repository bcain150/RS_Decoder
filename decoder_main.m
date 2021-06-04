%{ 
Brendan Cain
Error Correcting Codes Proj 4
Reed-Solomon Code Decoder

The main MATLab script for testing Berlekamp-Massey Algorithm functions
used for decoding Reed-Solomon codeword in GF(2^m)
%}

num_tests = 0;
failed_tests = 0;
succeed_tests = 0;
fail_msgs = "";

%HW 8 tests
syms x;
prim = x^3 + x + 1;
prim_poly = sym2poly(prim);
m = 3;
t = 2;
p = 2;
n = p^m - 1;
k = n - 2*t;
disp("TEST 1 : Given the information below, complete two tests.");
disp("  Let p(x) be the primitive polynomial for the implementation of a narrow sense");
disp("  (7,3) Reed-Solomon code shown below:");
fprintf("  Polynomial form:\n\t");
print_poly("p(x)", prim_poly, false);
fprintf("  Power form:\n\t");
print_poly("p(x)", prim_poly, true);
fprintf("\n  Generating mapping between binary groups in GF(2^3) and powers of alpha for arithmetic...\n");
GF2m = GenerateGF2(m, prim_poly);

fprintf("  Given generator G(x) and message M(x) (both shown below), encode the codeword to send.\n");

G = [0 3 0 1 3]; % g(x) = x^4 + a^3x^3 + x^2 + ax + a^3
M = [0 4 2]; % m(x) = x^2 + a^4 + a^2

fprintf("  Polynomial forms:\n\t");
print_poly("G(x)", G, false);
fprintf("\t");
print_poly("M(x)", M, false);
fprintf("  Power forms:\n\t");
print_poly("G[x]", G, true);
fprintf("\t");
print_poly("M[x]", M, true);
fprintf("\n  Then after encoding (generator method), the codeword C(x) is created\n");
fprintf("  (Shown in test cases)\n");

% generator method of encoding
C = PolyMultGF2(G, M, GF2m);
R_err = [0 6 5 3 -1 5 5];

num_tests = num_tests+1;
try
    do_test("TEST 1.1", R_err, GF2m, G, C, M, true);
    succeed_tests = succeed_tests + 1;
catch test_err
    warning("%s", test_err.message);
    failed_tests = failed_tests + 1;
    fail_msgs = [fail_msgs;test_err.message];
end


fprintf("\n  Now complete the second test...");
R_err = [0 6 5 3 -1 5 1];

num_tests = num_tests+1;
try
    do_test("TEST 1.2", R_err, GF2m, G, C, M, true);
    succeed_tests = succeed_tests + 1;
catch test_err
    warning("%s", test_err.message);
    failed_tests = failed_tests + 1;
    fail_msgs = [fail_msgs;test_err.message];
end

%HW 7 based test
fprintf("\nTEST 2 : Given the information below, complete three tests.\n");
prim = x^4 + x^3 + 1;
p = 2;
m = 4;
t = 3; %three error correcting code
n = p^m - 1; %codeword length
k = n - 2*t; %message length
prim_poly = sym2poly(prim); %array representation of prim variable

disp("  Let p(x) be the primitive polynomial for the implementation of a narrow sense");
disp("  (15,9) three-error correcting Reed-Solomon code shown below:");
fprintf("  Polynomial form:\n\t");
print_poly("p(x)", prim_poly, false);
fprintf("  Power form:\n\t");
print_poly("p[x]", prim_poly, true);
fprintf("\n  Generating mapping between binary groups in GF(2^4) and powers of alpha for arithmetic...\n");
GF2m = GenerateGF2(m, prim_poly);
fprintf("  Constructing the generator polynomial G(x)...\n");

%creating the generating polynomial using the second method
%the factors of G to be multiplied together in power form
%the placements of elements are powers of x, the value of the elements are
%the powers of alpha that are the x's coefficients
G = {[0 1], [0 2], [0 3], [0 4], [0 5], [0 6]};
while(size(G,2) > 1)
    G = {PolyMultGF2(G{1}, G{2}, GF2m), G{3:end}}; %appends G matrix
end
%G is now the generating polynomial for the RS code (Made using the 2nd method)
G = [G{:}];
fprintf("  Polynomial form:\n\t");
print_poly("G(x)", G, false);
fprintf("  Power form:\n\t");
print_poly("G[x]", G, true);

fprintf("\n  The next 2 tests utilize encoding the codeword via the generating\n");
fprintf("  method and the systematic method. The message polynomial is shown below:\n");

M = [1 -1 -1 -1 5 -1 -1 14 0]; %message polynomial

fprintf("\n  ----Generating Method----\n");
fprintf("  The generating the codeword to send (shown in test case)...\n");

C = PolyMultGF2(G, M, GF2m);

fprintf("  Constructing the simulated received codeword with 3 errors...\n");
R_err = [1 13 1 3 3 7 14 5 6 2 8 5 -1 13 6];

num_tests = num_tests + 1;
try
    do_test("TEST 2.1", R_err, GF2m, G, C, M, true);
    succeed_tests = succeed_tests + 1;
catch test_err
    warning("%s", test_err.message);
    failed_tests = failed_tests + 1;
    fail_msgs = [fail_msgs;test_err.message];
end

fprintf("\n  ----Systematic Encoding Method----\n");
shift = [0 -1 -1 -1 -1 -1 -1]; %this equals x^6
shifted = PolyMultGF2(M, shift, GF2m); %this is the shifted value after the computation
[~, rem] = PolyDivGF2(shifted, G, GF2m);
C_sys = shifted;

fprintf("  The systematically encoding the codeword to send (shown in test case)...\n");

C_sys(1,end-(n-k)+1:end) = rem; %systematic method codeword

fprintf("  Constructing the simulated received codeword with 3 errors...\n");
R_err = [1 -1 0 -1 5 2 -1 14 0 2 7 -1 5 3 13 ];

num_tests = num_tests + 1;
try
    do_test("TEST 2.2", R_err, GF2m, G, C_sys, M, false);
    succeed_tests = succeed_tests + 1;
catch test_err
    warning("%s", test_err.message);
    failed_tests = failed_tests + 1;
    fail_msgs = [fail_msgs;test_err.message];
end

fprintf("\n  ----Too Many Errors----\n");
fprintf("  Here, all variables are the same as the last test case except\n");
fprintf("  There are 4 errors in R(x) instead of 3. The decoder should fail.\n");
fprintf("  Constructing the simulated received codeword with 4 errors...\n");
R_err = [1 -1 0 -1 5 2 -1 14 0 2 7 -1 5 3 0];

num_tests = num_tests + 1;
try
    do_test("TEST 2.3", R_err, GF2m, G, C_sys, M, false);
    succeed_tests = succeed_tests + 1;
catch test_err
   warning("%s", test_err.message);
   failed_tests = failed_tests + 1;
   fail_msgs = [fail_msgs;test_err.message];
end

fprintf("\n\n-------------TESTING COMPLETE-------------\n");
fprintf("Statistics:\n");
fprintf("\tNumber of Tests: %d\n", num_tests);
fprintf("\tSucceeded: %d\n", succeed_tests);
fprintf("\tFailed: %d\n", failed_tests);
fprintf("Failure Messages:\n");

if(size(fail_msgs,1) < 2)
    fprintf("\tThere are no messages to display\n");
else
    for i=2:size(fail_msgs,1)
        fprintf("%d.) %s\n", i-1, fail_msgs(i));
    end
end



