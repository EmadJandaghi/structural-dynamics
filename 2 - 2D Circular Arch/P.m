function load=P(t)
    global nDOFs
    P1=0.6*(1e6); % N
    td=0.05;
    load=zeros(nDOFs,1);
    if t<=td
        load(2*32)=(P1/td)*t;
    else
        load(2*32)=P1;
    end
    
    
    
end
