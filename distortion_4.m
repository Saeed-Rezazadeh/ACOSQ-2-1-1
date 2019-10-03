function Distortion = distortion_3(f ,  y_1 , y_2 , y_3 , codebook , delta , Pr_z , T)
y_1_2 = (y_1 - 1) * 2 + y_2 ;
summation = 0 ;
parfor x_4 = 1 : 2
    for y_4 = 1 : 2
        u_index = find (T(: , 7) == x_4) ;
        for u_i = 1 : length(u_index)
            x_3 = T(u_index(u_i) , 2 + y_1_2) ; 
            summation = summation + Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1 )...
                * delta * f(u_index(u_i)).* (T(u_index(u_i) , 1) - codebook(y_4)) .^ 2 ;
        end
    end
end
Distortion = summation ;
end