eseS = dot(es,eS,4); 
eueU = dot(eu,eU,4);
eveV = dot(ev,eV,4);

eseS_is_one = ismembertol(eseS,ones(size(eseS)));
eueU_is_one = ismembertol(eueU,ones(size(eueU)));
eveV_is_one = ismembertol(eveV,ones(size(eveV)));

%% plot eseS

for i=1:data.ns
    d1 = reshape(eseS_is_one(i,:,:),[size(u),size(v)]);
    if (min(d1) == 0 && i ~= 1)
        figure()
        pcolor(eseS(i,~d1))
        xlabel('u')
        ylabel('v')
        title(sprintf('es $\cdot$ eS at s =%f',data.phi(i)))
    end
end

%% plot eueU

for i=1:data.ns
    d1 = reshape(eueU_is_one(i,:,:),[size(u),size(v)]);
    if (min(d1) == 0 && i ~= 1)
        figure()
        pcolor(eueU(i,~d1))
        xlabel('u')
        ylabel('v')
        title(sprintf('eu $\cdot$ eU at s =%f',data.phi(i)))
    end
end

%% plot eveV

for i=1:data.ns
    d1 = reshape(eveV_is_one(i,:,:),[size(u),size(v)]);
    if (min(d1) == 0 && i ~= 1)
        figure()
        pcolor(eveV(i,~d1))
        xlabel('u')
        ylabel('v')
        title(sprintf('ev $\cdot$ eV at s =%f',data.phi(i)))
    end
end

% check cross terms
