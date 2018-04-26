% Written by Marco Pingaro and Paolo Venini
%
function er_u = error_l2_norm_rot(sp, mc, el, cr, l)
    nel = size(el,1);
    
    [gauss_w, gauss_p] = GaussQuad2D(4,4);
    npg = size(gauss_w,1);
   
    err_u = 0.;
    %norm_u = 0.;
    
    for i=1:nel
        pt = [cr(el(i,1),1),cr(el(i,1),2);
              cr(el(i,2),1), cr(el(i,2),2);
              cr(el(i,3),1), cr(el(i,3),2);
              cr(el(i,4),1), cr(el(i,4),2)];
        
        for j=1:npg
            
            [DF,DFF,JF] = jacobian_quad(pt,gauss_p(j,1),gauss_p(j,2));
            [x,y] = map_quad(pt,gauss_p(j,:));

            %sol = 0.5*pi*cos(pi*x)*( 2*cos(2*pi*x) - cos(pi*y) );
            %sol = (pi*cos(pi*x)*( cos(pi*y) - 4*cos(pi*y)^2 + 2))/2;
            sol =-x^2*(x-1)^2*(2*y^2-3*y+1)-y^2*(y-1)^2*(2*x^2-3*x+1)-x*y^2*(4*x-3)*(y-1)^2-x^2*y*(4*y-3)*(x-1)^2;

            rot(1,1) = 0.25*(1-gauss_p(j,1))*(1-gauss_p(j,2)) ;
            rot(1,2) = 0.25*(1+gauss_p(j,1))*(1-gauss_p(j,2)) ;
            rot(1,3) = 0.25*(1+gauss_p(j,1))*(1+gauss_p(j,2)) ;
            rot(1,4) = 0.25*(1-gauss_p(j,1))*(1+gauss_p(j,2)) ;                         
            
            sp_r =  sp(el(i,1),1)*rot(1,1) + sp(el(i,2),1)*rot(1,2) +...
                sp(el(i,3),1)*rot(1,3) + sp(el(i,4),1)*rot(1,4);
            
            err_u  = err_u + ((sp_r - sol)^2)*gauss_w(j)*JF;
            %norm_u = norm_u + (sol^2)*gauss_w(j);

        end
    end
    %
    %npt = nel*npg;
    %norm_u = sqrt ( norm_u ) ;
    er_u = sqrt ( err_u );
    
end
