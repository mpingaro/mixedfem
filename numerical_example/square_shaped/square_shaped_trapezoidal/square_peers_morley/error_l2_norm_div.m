% Written by Marco Pingaro and Paolo Venini
%
function er_u = error_l2_norm_div(sp, mc, el, cr, l, mu)
    nel = size(el,1);
    
    [gauss_w, gauss_p] = GaussQuad2D(4,4);
    npg = size(gauss_w,1);
   
    err_ux = 0.;
    err_uy = 0.;
    %norm_ux = 0.;
    %norm_uy = 0.;
    
    for i=1:nel
        pt = [cr(el(i,1),1),cr(el(i,1),2);
              cr(el(i,2),1), cr(el(i,2),2);
              cr(el(i,3),1), cr(el(i,3),2);
              cr(el(i,4),1), cr(el(i,4),2)];
        
        s = mc(i,13) ;
        
        for j=1:npg
            
            xi  = gauss_p(j,1);
            eta = gauss_p(j,2);
            [DF,DFF,JF] = jacobian_quad(pt,xi,eta);
            [x,y] = map_quad(pt,gauss_p(j,:));
            
            %sol(1,1) = -pi^2*cos(pi*x)*sin(pi*y)*(l + mu + 2*l*cos(pi*y) +...
            %    12*mu*cos(pi*y));
            %sol(1,2) = -pi^2*sin(pi*x)*(l*cos(pi*y) + 3*mu*cos(pi*y) +...
            %    2*l*(2*cos(pi*y)^2 - 1) + 2*mu*(2*cos(pi*y)^2 - 1));

            sol(1,1) = -4*mu*(2*y-1)*(3*x^4-6*x^3+6*x^2*y^2-6*x^2*y+3*x^2-6*x*y^2+6*x*y+y^2-y);
            sol(1,2) = 4*mu*(2*x-1)*(6*x^2*y^2-6*x^2*y+x^2-6*x*y^2+6*x*y-x+3*y^4-6*y^3+3*y^2);

            div(:,1)  = s*[0.5, 0] ; 
            div(:,2)  = s*[0, 0.5] ; 
            div(:,3)  = s*[0.5, 0] ; 
            div(:,4)  = s*[0, 0.5] ; 
            div(:,5)  = s*[0.5, 0] ; 
            div(:,6)  = s*[0, 0.5] ; 
            div(:,7)  = s*[0.5, 0] ; 
            div(:,8)  = s*[0, 0.5] ; 
            div(:,9) = [0, 0] ; 
            div(:,10) = [0, 0] ;
            div(:,11) = [0, 0] ;
            div(:,12) = [0, 0] ;
            
            sp_divx = sp(mc(i,1:12))'*div(1,:)';
            sp_divy = sp(mc(i,1:12))'*div(2,:)';
            
            err_ux  = err_ux + ((sp_divx/JF - sol(1,1))^2)*gauss_w(j)*JF;
            err_uy  = err_uy + ((sp_divy/JF - sol(1,2))^2)*gauss_w(j)*JF;
            
            %norm_ux = norm_ux + (sol(1,1)^2)*gauss_w(j);
            %norm_uy = norm_uy + (sol(1,2)^2)*gauss_w(j);
            
        end
    end
    %
    %npt = nel*npg;
    %norm_u = sqrt ( norm_ux + norm_uy ) ;
    er_u = sqrt ( (err_ux + err_uy) ); %/norm_u;
   
end
