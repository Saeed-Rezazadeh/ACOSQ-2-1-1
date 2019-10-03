function [SDR , Distortion , T , codebook] = COSQ_3(f , y_1 , y_2 , Pr_z , T , codebook , delta )

FileID = fopen ('Results.txt' , 'a') ;

D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal Partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        temp = zeros(2 , 1) ;
        u = T(u_index , 1) ;
        x_prime = T(u_index , 2) ;
        if (mod(x_prime - 1 , 2) == 0)
            x_2 = 1 ;
        else
            x_2 = 2 ;
        end
        for x_3 = 1 : 2
            for y_3 = 1 : 2
                summation = summation + Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) * (u - codebook(y_3)) ^ 2 ;
            end
            temp (x_3) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min(temp) ;
        T_u(u_index , 1) = partition_index ;
    end
    T (: , 3) = T_u ;
    
    %% Optimal Centroids
    parfor y_3 = 1 : 2
        numerator = 0 ;
        denominator = 0 ;
        for x_3 = 1 : 2
            u_index = find (T(: , 3) == x_3) ;
            for u_i = 1 : length(u_index)
                x_prime = T(u_index(u_i) , 2) ;
                if (mod(x_prime - 1 , 2) == 0)
                    x_2 = 1 ;
                else
                    x_2 = 2 ;
                end
                
                numerator = numerator + ... 
                    Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ... 
                    * T(u_index(u_i) , 1) * f(u_index(u_i)) ;
                denominator = denominator + ... 
                    Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1 ) ...
                    * f(u_index(u_i)) ;
            end
        end
        codebook(y_3) = numerator / denominator ;
    end
    %% Distortion
    D(2) = distortion_3(f ,  y_1 , y_2 , codebook , delta , Pr_z , T) ;
    fprintf (FileID , 'Overall D_3 = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_3 = %f\n' , SDR) ;
fclose (FileID) ;
end