%/***************************************************************************************/
function Jtime = tmtonumb_d( yr, xmth, day, hr, min, sec )
%/***************************************************************************************/
 days(1:12)=[31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.];

jtime(1:length(yr)) = 0.0;
xjday(1:length(yr)) = 0.0;
jtime = jtime';
xjday = xjday';
for j=1:length(yr),
    xjday(j) = sum(days(1:xmth(j)-1)',1);
end
xjday = xjday + day;
if( mod( yr , 4. ) == 0 )
    if(xmth > 2) 
        xjday= xjday+1.; 
    end
end
Jtime = xjday + ( hr + min /60. + sec / 3600. ) / 24.;
return
