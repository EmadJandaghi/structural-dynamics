function load=P(t)
    global nDOFs
    loadedNodes=[61:90,121:150,181:210,241:270,301:330,361:390];
    loadedDOFs=3*loadedNodes;
    P1=-250; % N
    wb=250;
    load=zeros(nDOFs,1);
    load(loadedDOFs)=P1*sin(wb*t);

    
    
    
end
