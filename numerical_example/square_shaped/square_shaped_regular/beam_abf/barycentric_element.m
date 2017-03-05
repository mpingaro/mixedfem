function [B, A] = barycentric_element(Element, Node)

    nelem = size(Element,1);
    
    B = zeros(nelem,2);
    A = zeros(nelem,1);

    for i=1:nelem
        nl = numel(Element(i,:));
        
        % Computed Area of elements
        p = zeros(nl,2);
        for k=1:nl
            p(k,[1 2]) = [Node(Element(i,k),1),Node(Element(i,k),2)];
        end
        A(i) = area_polygon(p,nl);
        
        for j=1:nl
            
            if (j==nl)
                
                B(i,1) = B(i,1) + 0.5*( Node(Element(i,1),2) - Node(Element(i,j),2) ) * ...
                    ( (( Node(Element(i,1),1) - Node(Element(i,j),1 ))^2)/3 + Node(Element(i,j),1)*Node(Element(i,1),1) )/A(i);
                B(i,2) = B(i,2) - 0.5*( Node(Element(i,1),1) - Node(Element(i,j),1) ) * ... 
                    ( (( Node(Element(i,1),2) - Node(Element(i,j),2 ))^2)/3 + Node(Element(i,j),2)*Node(Element(i,1),2) )/A(i);
       
            else
     
                B(i,1) = B(i,1) + 0.5*( Node(Element(i,j+1),2) - Node(Element(i,j),2) ) * ... 
                    ( (( Node(Element(i,j+1),1) - Node(Element(i,j),1 ))^2)/3 + Node(Element(i,j),1)*Node(Element(i,j+1),1) )/A(i);
                
                B(i,2) = B(i,2) - 0.5*( Node(Element(i,j+1),1) - Node(Element(i,j),1) ) * ...
                    ( (( Node(Element(i,j+1),2) - Node(Element(i,j),2 ))^2)/3 + Node(Element(i,j),2)*Node(Element(i,j+1),2) )/A(i);
            end
        end
    end
    
