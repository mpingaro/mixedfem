% Written by Marco Pingaro and Paolo Venini
%
function er_u = error_l2_norm(sp, mc, el, cr, l)
    nel = size(el,1);
    sp = reshape(sp, 2, [])';
    
    [gauss_w, gauss_p] = GaussQuad2D(4,4);
    npg = size(gauss_w,1);
   
    err_ux = 0.;
    err_uy = 0.;
    norm_ux = 0.;
    norm_uy = 0.;
    
    for i=1:nel
        pt = [cr(el(i,1),1),cr(el(i,1),2);
              cr(el(i,2),1), cr(el(i,2),2);
              cr(el(i,3),1), cr(el(i,3),2);
              cr(el(i,4),1), cr(el(i,4),2)];
        
        for j=1:npg
            
            [DF,DFF,JF] = jacobian_quad(pt,gauss_p(j,1),gauss_p(j,2));
            [x,y] = map_quad(pt,gauss_p(j,:));
            
            sol_x = cos(pi*x)*sin(2*pi*y);
            sol_y = sin(pi*x)*cos(pi*y);
                         
            sp_x = sp(i,1);
            sp_y = sp(i,2);
             
            err_ux  = err_ux + ((sp_x - sol_x)^2)*gauss_w(j)*JF;
            err_uy  = err_uy + ((sp_y - sol_y)^2)*gauss_w(j)*JF;
            
            norm_ux = norm_ux + (sol_x^2)*gauss_w(j);
            norm_uy = norm_uy + (sol_y^2)*gauss_w(j);

        end
    end
    %
    norm_u = sqrt ( norm_ux + norm_uy ) ;
    er_u = sqrt ( err_ux + err_uy )/norm_u;
    
end