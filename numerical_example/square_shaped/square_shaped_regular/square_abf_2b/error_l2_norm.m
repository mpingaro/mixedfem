% Written by Marco Pingaro and Paolo Venini
%
function er_u = error_l2_norm(sp, el, cr, l)
    nel = size(el,1);

    [gauss_w, gauss_p] = GaussQuad2D(8,8);
    npg = size(gauss_w,1);
    
    sp = reshape(sp, 6, []);    
    npt = nel*npg;
    
    sol_a = zeros(npt,2);
    sp_c  = zeros(npt,2);
    %[bari, area] = barycentric_element(el, cr);


    for i=1:nel
        for j=1:npg
            pt = [cr(el(i,1),1),cr(el(i,1),2);
                cr(el(i,2),1), cr(el(i,2),2);
                cr(el(i,3),1), cr(el(i,3),2);
                cr(el(i,4),1), cr(el(i,4),2)];
        
            [x,y] = map_quad(pt,gauss_p(j,:));
            
            %A = 2/(1+l);
            %B = 0.5*A*sin(pi*pt(1))*sin(pi*pt(2));
            %b = 1/25;

            %sol_a(i,1) = b*(sin(2*pi*pt(2))*(-1+cos(2*pi*pt(1))) + B);
            %sol_a(i,2) = b*(sin(2*pi*pt(1))*( 1-cos(2*pi*pt(2))) + B);
        
            sp_c(npg*(i-1)+j,1)  = sp(1,i) + sp(2,i)*gauss_p(j,1) + sp(3,i)*gauss_p(j,2);
            sp_c(npg*(i-1)+j,2)  = sp(4,i) + sp(5,i)*gauss_p(j,1) + sp(6,i)*gauss_p(j,2);

            sol_a(npg*(i-1)+j,1) = -2*x^2*y*(1-x)^2*(1-3*y+2*y^2);
            sol_a(npg*(i-1)+j,2) =  2*x*y^2*(1-y)^2*(1-3*x+2*x^2);
        end
    end
    %er_ux = sum( (sp(:,1) - sol_a(:,1) ).^2 )/nel ;
    %er_uy = sum( (sp(:,2) - sol_a(:,2) ).^2 )/nel ;
    
    er_ux = sum( (sp_c(:,1) - sol_a(:,1) ).^2 )/npt ;
    er_uy = sum( (sp_c(:,2) - sol_a(:,2) ).^2 )/npt ; 
    %
    norm_u = sqrt ( sum( sol_a(:,1).^2 + sol_a(:,2).^2 ) );
    er_u = sqrt ( er_ux + er_uy )/norm_u;
end            
