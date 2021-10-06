# Perterbed Lambert problem solver        
# Reference: Engels & Junkins, The gravity-perturbed Lambert probelm: A KS Variation of Parameters Approach, (1980)

function lambert(R0,Rf,T;
    R = 6378.137, # Earth radius (km)
    mu = 398600.8, # Earth graviational constant (km^3/sec^2)
    J2 = 0.0010826157) #set to 0 for unperturbed
    
    # Inputs: 
    #   R0 (float64,3-vector) -> initial position vector (km)
    #   Rf (float64,3-vector) -> final position vector (km)
    #   T (float64) -> time of flight (sec)
    #   Optional:
    #       R (float64) -> equatorial radius of gravitational body
    #       mu (float64) -> graviational constant of gravitational body
    #       J2 (float64) -> value of the J2 constant

    
    #(1) Given: r(0),r(f),t(f)-t(0)
    # Normalize to canonical units    
    km2cdu = 1/R
    R0 = R0*km2cdu/1000 #cdu
    Rf = Rf*km2cdu/1000 #cdu
    
    s2ctu = 1/sqrt(R^3/mu) #sec
    T = T*s2ctu #ctu

    print("R0 = $R0, Rf = $Rf, T = $T\n")
    #(2) Construct u(0) from equations (A2-A3) and u0(sf0) from equations (A12-A13)   
    # construct u(0)
    r0 = sqrt(R0'*R0)
    x0,y0,z0 = R0
    u0 = zeros(4)
    if x0>=0#(A2)
        u0[1] = sqrt((r0+x0)/2)
        u0[2] = (y0*u0[1])/(r0+x0)
        u0[3] = (z0*u0[1])/(r0+x0)
        u0[4] = 0
    else #(A3)
        u0[2] = sqrt((r0-x0)/2)
        u0[1] = (y0*u0[2])/(r0-x0)
        u0[3] = 0
        u0[4] = (z0*u0[2])/(r0-x0)
    end

    # construct u0(sf0), denoted by uf0 in code
    rf = sqrt(Rf'*Rf)
    xf,yf,zf = Rf

    uf0 = zeros(4)
    if xf >= 0 #(A12)
        P = (u0[4]*(rf+xf)+u0[2]*zf-u0[3]*yf)/(u0[1]*(rf+xf)+u0[2]*yf+u0[3]*zf)

        uf0[1] = sqrt((rf+xf)/2*(1+P^2))
        uf0[4] = P*uf0[1]
        uf0[2] = (yf*uf0[1]+zf*uf0[4])/(rf+xf)
        uf0[3] = (zf*uf0[1]-yf*uf0[4])/(rf+xf)
    else # (A13)
        Q = (u0[2]*(rf-xf)+u0[1]*yf+u0[4]*zf)/(u0[3]*(rf-xf)+u0[1]*zf-u0[4]*yf)

        uf0[3] = sqrt((rf-xf)/(2*(1+Q^2)))
        uf0[2] = Q*uf0[3]
        uf0[1] = (zf*uf0[3]+yf*uf0[2])/(rf-xf)
        uf0[4] = (zf*uf0[2]-yf*uf0[3])/(rf-xf)
    end
    
    #(3) Solve the unperturbed universal Lambert problem using Equations (28-37) with u(f)= u0(sf0), yielding at0,sf0, and udot0(0). This involves an iterative process.    
    t0 = 0
    tf = T

    #at0,sf0 = atsf(r0,rf,u0,uf0,t0,tf)

    epsilon = 1f-8
    ctr = 0    
    iterationlimit = 1000    
    at0 = 1 #first guess
    sf0 = 1  #first guess
    del = [1;1] #kickstart

    while any(abs.(del).> epsilon)
        ctr = ctr+1        
        if ctr >= iterationlimit
            print("Could not converge in $ctr iterations\n")
            break                        
        end        

        # update stumpff functions for shorthand notation
        c(n) = stumpff(n,at0*sf0^2)
        ct(n) = stumpff(n,4*at0*sf0^2)
        
        # F,G 2 eqations in 2 uknowns at,sf
        F = rf + r0 - sf0^2 * ct(2) - 2 * u0' * uf0 * c(0) #(30)
        G = tf - t0 - sf0^3 * ct(3) - u0' * uf0 * sf0 * c(1) #(31)  

        # partials for Jacobian in Newton's method
        Fpsf = sf0 * c(1) * (2 * at0 * u0' * uf0 - c(0)) #(32)
        Gpsf = -sf0^2 * ct(2) - u0' * uf0 * c(0) #(33)
        Fpat = sf0^2 * (u0' * uf0 * c(1) - 2 * sf0^2 * (2 * ct(4) - ct(3))) #(34)
        Gpat = -sf0^3 * (2 * sf0^2 * (3 * ct(5) - ct(4)) + (1/2) * u0' * uf0 * (c(3) - c(2))) #(35)
        D = (Fpat*Gpsf-Fpsf*Gpat)
        
        # Newton's method
        tmp = [at0;sf0]        
        (at0,sf0) = [at0;sf0] - (1/D)*[Gpsf -Fpsf; -Gpat Fpat]*[F;G]
        print("at0 = $at0, sf0 = $sf0\n")
        del = ([at0;sf0] - tmp)
    end    
    if ctr < iterationlimit
        print("Converged in $ctr iterations\n")        
    end   
    
    # Solve for udot0 and rdot0 (v0)
    cc(n) = stumpff(n,at0*sf0^2)    #cc to prevent warning for overriding previous method
    udot0 = 1/(sf0*cc(1))*(uf0-u0*cc(0)) #(37)
    B = 2/r0*[u0[1] -u0[2] -u0[3] u0[4]; u0[2] u0[1] -u0[4] -u0[3]; u0[3] u0[4] u0[1] u0[2]] #(39)
    v0 = B*udot0 #(38)    
    v0 = v0*(s2ctu/km2cdu) #convert back from canonical units to metric units

 #=
    #Solve the integrals in equations 76-79
    V1 = 1/2((3*zf^2-rf^2)/rf^5) #(B1)

    M = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]
    I = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    K = 1/rf^5*(3*rf*zf*M-(6*zf^2-rf^2)*I)
    Q1 = -1/2*K*u0
    G1 = u0'*Q1

    alpha1 = quadgk(s -> Q1*s*stumpff(1,at0*s^2), 0, sf0, rtol=1e-8)
    beta1 = quadgk(s -> Q1*stumpff(0,at0*s^2), 0, sf0, rtol=1e-8)

    # Solve for q1
    l(v,w) = v[4]*w[1]-v[3]*w[2]+v[2]*w[3]-v[1]*w[4]    
    q1 = l(u0, (-alpha1*stumpff(0,at0*sf0^2)+beta1*sf0*stumpff(1,at0*sf0^2))) #(A18)

    #Solve for u1 (A21)
    P = 2*(u0[1]^2+u0[4]^2)*(u0[4]*uf0[4]+u0[1]*uf0[1])+yf*(u0[2]*uf0[1]-u0[3]*uf0[4])+zf*(u0[2]*uf0[4]+u0[3]*uf0[1])

    u1[1] = 2*q1*(u0[1]^2+u0[4]^2)*u0[4]/P
    u1[4] = -uf0[1]/uf0[4]*u1[1]
    u1[2] = (yf*u1[1]+zf*u1[4])/(2*(uf0[1]^2+uf0[4]^2))
    u1[3] = (zf*u1[1]-yf*u1[4])/(2*(uf0[1]^2+uf0[4]^2))

    #Solve for tau1 and a1

    #partials
    pps_c0 = 0
    ppat_c0 = 0
    pps_sc1 = stumpff(0,0)
    ppat_sc1 = 0
    pps_s2ct2 = stumpff(0,0)
    ppat_s2ct2 = stumpff(0,0)
    pps_s3ct3 = stumpff(0,0)
    ppat_s3ct3 = stumpff(0,0)
 =#

    return v0
end

function stumpff(k,x) 
    # calculates Stumpff functions c, ct 
    # Karl Stumpff (1956), Victor Bond (1974)
    if x == 0
        c = 1/factorial(k)
    else
        tmp = zeros(k+1)
        for i = 0:k
            if i == 0 
                if x > 0
                    tmp[i+1] = cos(sqrt(x))
                elseif x < 0
                    tmp[i+1] = cosh(sqrt(-x))                    
                end                
            elseif i == 1
                if x > 0
                    tmp[i+1] = sin(sqrt(x))/sqrt(x)
                elseif x < 0 
                    tmp[i+1] = sinh(sqrt(-x))/sqrt(-x)                    
                end                
            else
                tmp[i+1] = (1/factorial(i-2) - tmp[i+1-2])/x                
            end
        end
        c = tmp[end]
    end
return c
end