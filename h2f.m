function f = h2f(var,ns)
% Map quantitiy from half to full grid
if size(var) < ns
    temp=zeros(1,ns);
    temp(1)      =  1.5 *  var( 1) - 0.5 *    var(2);
    temp(2:ns-1) =  0.5 * ( var(1:ns-2) + var(2:ns-1));
    temp(ns)     =  1.5*      var(ns-1) - 0.5* var(ns-2);
    f=temp;

else
    temp=zeros(1,ns);
    temp(1)      =  1.5 *   var( 2) - 0.5 *    var(3);
    temp(2:ns-1) =  0.5 * ( var(2:ns-1) + var(3:ns));
    temp(ns)     =  2   *   var(ns) - var(ns-2);
    f=temp;
end
end
