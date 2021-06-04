function print_poly(name, polynomial, pwr_form)
%{
PRINT_POLY: prints the given polynomial to the command window in the form
            requested
INPUTS:
    name - a string containing the name of the polynomial ex: "f(x)"
    polynomial - the power form representation of that polynomial
    pwr_form - a boolean variable that tells the function what verison to
               print. If true, power form representation printed. If false,
               polynomial form is printed
%}


poly_len = size(polynomial,2);

if(pwr_form)
    print_string = sprintf("%s = [", name);
else
    print_string = sprintf("%s = ", name);
end

for i = 1:poly_len
    curr_a = polynomial(i);
    %if power form
    if(curr_a == -1 && pwr_form)
        print_string = strcat(print_string, " Inf");
    elseif(pwr_form)
        curr_string = sprintf(" %d", curr_a);
        print_string = strcat(print_string, curr_string);
    %if polynomial form
    elseif(curr_a ~= -1)
        x_pwr = poly_len - i;
        if(curr_a == 1)
            a_str = "a";
        else
            a_str = sprintf("a^%d", curr_a);
        end
        if(x_pwr == 1)
            x_str = "x";
        else
            x_str = sprintf("x^%d", x_pwr);
        end
        
        if(curr_a ~= 0 && x_pwr ~= 0)
            curr_string = sprintf("%s*%s + ", a_str, x_str);
        elseif(x_pwr ~= 0)
            curr_string = sprintf("%s + ", x_str);
        elseif(curr_a ~= 0)
            curr_string = sprintf("%s + ", a_str);
        else
            curr_string = "1 + ";
        end
        print_string = strcat(print_string, curr_string);
    end    
end

if(pwr_form)
    fprintf("%s ]\n", print_string);
else
    len = strlength(print_string);
    print_string = eraseBetween(print_string,len-2, len);
    fprintf("%s\n", print_string);
end


end

