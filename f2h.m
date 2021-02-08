function f = f2h(var,ns)
% Map quantitiy from full to half grid
s = size(var);
temp=zeros(s(1),ns-1);
temp(:,:)      =  0.5 *   (var(:, 2:ns ) +    var(:,1:ns-1));
f=temp;
end
