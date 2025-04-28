function [ind_min_max]=get_min_max(inarr, Tolerance)
%==================================
Npoints=length(inarr);
k=2;
while k<=Npoints
    if ( abs(inarr(k)-inarr(1)) > Tolerance )
        break
    end
    k=k+1;
end
if( k > Npoints)
    k=Npoints;
end
ln_ind(1) = 1;
minmax=2;
if( inarr(k) > inarr(1) )
    trace = 'MAX';
    maxz = inarr(k);
    ln_ind(minmax) = k;
else
    trace = 'MIN';
    minz = inarr(k);
    ln_ind(minmax) = k;
end

for j=k:Npoints
    if( strcmp(trace, 'MAX') > 0 )
        if( inarr(j) >= maxz )
            maxz = inarr(j);
            ln_ind(minmax) = j;
        else
            if( Tolerance <= abs(maxz - inarr(j)) )
                minz = inarr(j);
                minmax=minmax+1;
                ln_ind(minmax) = j;
                trace = 'MIN';
            end
        end
    else
        if( inarr(j) <= minz )
            minz = inarr(j);
            ln_ind(minmax) = j;
        else
            if( Tolerance <= abs(inarr(j) - minz) )
                maxz = inarr(j);
                minmax=minmax+1;
                ln_ind(minmax) = j;
                trace = 'MAX';
            end
        end
    end
end
%==================================

np=length(ln_ind);
if(ln_ind(np)<Npoints)
    ln_ind(np)=Npoints;
end

ind_min_max=ln_ind;

return
