% by Marco Pingaro

function def = defomesh(sp,elem,coord)

sp = sp( 1:3:size(sp,1)-2 ) ;
sp = reshape(sp,2,[])' ;

def = zeros( size(coord) ) ;
for i = 1:size(coord,1)
    
    ind = [];
    for j =1:size(elem,1)
        if ~isempty( find( elem(j,:) == i) )
            in = j ;
        else
            in = [] ;
        end
        ind = union(ind,in) ;
    end
    def(i,1) = coord(i,1) + sum( sp(ind,1) )/length(ind) ;
    def(i,2) = coord(i,2) + sum( sp(ind,2) )/length(ind) ;
end

return