using GalacticOptim, Optim

function objlambert(v,p)  
    r1 = p[1:3]  
    r2 = p[4:6]
    t1 = p[7]
    t2 = p[8]

    x = [r1;v]
    T = t2-t1
    sol = sc(x,T)
    rf = sol.u[end,1][1:3]

    er = rf - r2
    e = er'*er
    return e
end

function optlambert(D1,D2,t1,t2)    
    r1 = D1.r
    v0 = D1.v
    r2 = D2.r
    p = [r1;r2;t1*86400;t2*86400]   

    fobj = OptimizationFunction(objlambert, GalacticOptim.AutoForwardDiff())
    prob = OptimizationProblem(fobj,v0,p)
    sol = solve(prob,NewtonTrustRegion(),f_tol = 1e-9, x_tol = 1e-9)

    return sol
end

function dolambert(N1,N2,T1,T2)
        
    R1 = oe2rv(T1,N1)
    R2 = oe2rv(T2,N2)

    sol = optlambert(R1,R2,T1,T2)    

    x0 = [R1.r;sol.u]   
    tspan = (T1*86400,T2*86400) 
    prop = sc(x0,tspan)
    rf = prop.u[end][1:3]

    e = rf-R2.r
    emag = sqrt(e'*e)
    
    
    #calc delta v
    dv1 = sol.u - R1.v
    dv2 = R2.v - prop.u[end][4:6]    

    dv1mag = sqrt(dv1'*dv1)
    dv2mag = sqrt(dv2'*dv2)
    # calc mass change    
    ve = Isp*g0
    m1 = 5000*exp(-dv1mag/ve)
    m2 = 5000*exp(-dv1mag/ve)

    m = m1+m2

    return sol.u,emag,m#,prop
end






    