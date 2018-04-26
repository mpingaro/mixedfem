% Written by Marco Pingaro and Paolo Venini
%

function [er_ux er_uy] = error_l2_norm_div(sig, mc, el, cr, mu)
    nel = size(el,1);
    
    sol_a = seros(1,2)
    diff1 = 0.;
    diff2 = 0.;

    for i=1:size(el,1)
        for j=1:4 % Number of edges per element    
       
        if (j==4)
            x = ( cr( el(i,1),1 ) - cr( el(i,j),1 ) )/2.;        
            y = ( cr( el(i,1),2 ) - cr( el(i,j),2 ) )/2.;
        else
            x = ( cr( el(i,j+1),1 ) - cr( el(i,j),1 ) )/2.;        
            y = ( cr( el(i,j+1),2 ) - cr( el(i,j),2 ) )/2.;
        end


        sol_a(1) = -4*mu*(2*y - 1)*(3*x^4 - 6*x^3 + 6*x^2*y^2 - 6*x^2*y + 3*x^2 - 6*x*y^2 + 6*x*y + y^2 - y);
        sol_a(2) =  4*mu*(2*x - 1)*(6*x^2*y^2 - 6*x^2*y + x^2 - 6*x*y^2 + 6*x*y - x + 3*y^4 - 6*y^3 + 3*y^2); 

        diff1 = diff1 + (sol_a(1) - sig() )^2;
        diff2 = diff2 + (sol_a(2) - sig() )^2;
        end
    end
    er_div_1 = sqrt( sum( (sol_a(:,1) - sp(:,1)).^2 ) );
    er_div_2 = sqrt( sum( (sol_a(:,2) - sp(:,2)).^2 ) );
end            
