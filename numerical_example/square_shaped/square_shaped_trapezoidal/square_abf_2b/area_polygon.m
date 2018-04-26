function area = area_polygon(p,n)
    area = 0.0;
    j = n;       % The last vertex is the 'previous' one to the first

    for i=1:n
        area = area + ( p(j,1) + p(i,1) )*( p(j,2) - p(i,2) );
        j = i;  % j is previous vertex to i
    end
    area = abs(area/2.);