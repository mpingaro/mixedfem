% by Marco Pingaro

function def = defomesh(sp,elem,coord)

def = zeros( size(coord) ) ;
sp = reshape(sp,2,[])' ;

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
