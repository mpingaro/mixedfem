% Written by Marco Pingaro and Paolo Venini
%
function er_u = error_l2_norm(sp, el, cr, l)
    nel = size(el,1);

    sp = sp( 1:3:size(sp,1)-2 ) ;
    sp = reshape(sp, 2, [])';
    sol_a = zeros(nel,2);

    for i=1:size(el,1)
        pt(1) = ( cr(el(i,1),1) + cr(el(i,2),1) + cr(el(i,3),1) + cr(el(i,4),1) )/4;
        pt(2) = ( cr(el(i,1),2) + cr(el(i,2),2) + cr(el(i,3),2) + cr(el(i,4),2) )/4;

        A = 2/(1+l);
        B = 0.5*A*sin(pi*pt(1))*sin(pi*pt(2));
        b = 1/25;

        sol_a(i,1) = b*(sin(2*pi*pt(2))*(-1+cos(2*pi*pt(1))) + B);
        sol_a(i,2) = b*(sin(2*pi*pt(1))*( 1-cos(2*pi*pt(2))) + B);
    end
    er_ux = sum( (sp(:,1) - sol_a(:,1) ).^2 )/nel ;
    er_uy = sum( (sp(:,2) - sol_a(:,2) ).^2 )/nel ;
    %
    norm_u = sqrt ( sum( sol_a(:,1).^2 + sol_a(:,2).^2 ) );
    er_u = sqrt ( er_ux + er_uy )/norm_u;
end            
