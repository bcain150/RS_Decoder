function [bin, msg_symb] = get_message(code_wrd, GF, G, k)
% GET_MESSAGE Gets the message and bits of a codeword that was
%generated using the systematic or generator method. If the
%codeword cannot be decoded it throws an error
%INPUTS:
% code_wrd - the GF2^m representation of a message or codeword polynomial
% GF -  the array of cells enumerating GF(2^m) symbols as elements in GF(2)[a].
% G - represents the generating polynomial for a
%     code design so that code_wrd can be decoded (into binary) if its a
%     received polynomial.
% k - optional parameter, if k exists the encoding is systematic
%OUTPUTS:
% bin - the binary representation of code_wrd in GF(2) or the binary
%       representation of the message that code_wrd represents in GF(2)
% msg - the decoded message polynomial

gen_enc = false;
if(~exist('k', 'var'))
    gen_enc = true;
end

m = size(GF{1},2); %symbol length

if(gen_enc)
    [msg_symb, rem] = PolyDivGF2(code_wrd, G, GF);
    if(rem(:) ~= -1)
        error("Could not get the original message polynomial from the codeword");
    end
else
    msg_symb = code_wrd(1, 1:k);
end

wrd_len = size(msg_symb, 2);
binary_len = wrd_len*m;

for i = 1:wrd_len
    start = (i-1)*m;
    ending = i*m;
    if(ending > binary_len) 
        error("Binary index out of bounds.");
    end
    bin(1,1+start:ending) = GF{msg_symb(i)+2};
end



end

