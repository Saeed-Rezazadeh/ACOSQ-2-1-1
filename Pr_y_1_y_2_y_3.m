function Probability_y_1_2_3 = Pr_y_1_y_2_y_3( y_1, y_2 , y_3 , Pr , Pr_z , f  , T , delta)
y_1_2 = (y_1 - 1) * 2 + y_2 ;
summation = 0 ;
for x_1_2 = 1 : 4
    if (mod(x_1_2 - 1 , 2) == 0)
        x_2 = 1 ;
    else
        x_2 = 2 ;
    end
    for x_3 = 1 : 2
        u_index = find(T(: , 2) == x_1_2 & T(: , 2 + y_1_2) == x_3) ;
        
        summation = summation + Pr(x_1_2 , y_1_2) * Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
            * delta * sum(f(u_index));
        
    end
end
Probability_y_1_2_3 = summation ;
end