% Written by Marco Pingaro and Paolo Venini
%

function [er_ux er_uy] = error_l2_norm(sp, el, cr)
    nel = size(el,1);
    sp = reshape(sp, 2, [])';
    sol_a = zeros(nel,2);

    for i=1:size(el,1)
        pt(1) = ( cr(el(i,1),1) + cr(el(i,2),1) + cr(el(i,3),1) )/3;
        pt(2) = ( cr(el(i,1),2) + cr(el(i,2),2) + cr(el(i,3),2) )/3;

        %sol_a(i,1) = cos( pi*pt(1) )*sin( 2*pi*pt(2) );
        %sol_a(i,2) = sin( pi*pt(1) )*cos( pi*pt(2) );
        
        sol_a(i,1) = sin( pi*pt(1) )*sin( pi*pt(2) );
        sol_a(i,2) = sin( pi*pt(1) )*sin( pi*pt(2) );

    end
    er_ux = sqrt( sum( (sol_a(:,1) - sp(:,1)).^2 ) );
    er_uy = sqrt( sum( (sol_a(:,2) - sp(:,2)).^2 ) );
end            
