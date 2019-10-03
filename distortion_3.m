function Distortion = distortion_2(f ,  y_1 , y_2 , codebook , delta , Pr_z , T)
summation = 0 ;
parfor x_3 = 1 : 2
    for y_3 = 1 : 2
        u_index = find (T(: , 3) == x_3) ;
        for u_i = 1 : length(u_index)
            x_prime = T(u_index(u_i) , 2) ; 
            if (mod(x_prime - 1 , 2) == 0) 
                x_2 = 1 ; 
            else
                x_2 = 2 ; 
            end 
            summation = summation + Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 )...
                * delta * f(u_index(u_i)).* (T(u_index(u_i) , 1) - codebook(y_3)) .^ 2 ;
        end
    end
end
Distortion = summation ;
end