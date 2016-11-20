% Written by Marco Pingaro and Paolo Venini
%
% Solution present in the paper: Three field formulation

function er_u = error_beam_l2_norm(sp,el,cr,l,E,nu,f)
    np = size(el,1); % Number of element and evalution point
    
    sp = sp( 1:3:size(sp,1)-2 ) ;
    sp = reshape(sp, 2, [])';
    
    sol_a = zeros(np,2);

    for i=1:np
        pt(1) = (cr(el(i,1),1)+cr(el(i,2),1)+cr(el(i,3),1)+cr(el(i,4),1))/4 ;
        pt(2) = (cr(el(i,1),2)+cr(el(i,2),2)+cr(el(i,3),2)+cr(el(i,4),2))/4 ;


        sol_a(i,1) = 2*f*(1-nu^2)/(E*l)*pt(1)*(l/2 - pt(2));
        sol_a(i,2) = f/(E*l)*( pt(1)^2 + nu/(1-nu)*(pt(2)^2-l*pt(2)) );
    end
    er_ux = sum ( ( sp(:,1) - sol_a(:,1) ).^2 )/np;
    er_uy = sum ( ( sp(:,2) - sol_a(:,2) ).^2 )/np;

    norm_u = sqrt ( sum( sol_a(:,1).^2 + sol_a(:,2).^2 ) );
    er_u = sqrt ( er_ux + er_uy )/norm_u;
end 