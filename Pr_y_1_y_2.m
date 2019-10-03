function Probability_y_1_2 = Pr_y_1_y_2 (y_1 , y_2 , f , T , delta , Pr)
y = (y_1 - 1) * 2 + y_2 ;

summation = 0 ;
for x = 1 : 4
    u_index = find (T(: , 2) == x) ;
    summation = summation + Pr(x , y) * delta * sum(f(u_index)) ;
end
Probability_y_1_2 = summation ;
end