function [] = test(x)

aa = @(a) 2*a;

    function b = bb(a)
        b = 3*aa(a);
    end

disp( bb(x));

end