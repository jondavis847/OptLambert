using DataFrames, DifferentialEquations, LinearAlgebra

function oe2rv(t,o)
    #t  between   23467≤tevent≤26419
    oe = Data[o+1,:]
    t0,a,e,i,W0,w0,M0 = oe.t,oe.a,oe.e,oe.i,oe.W,oe.w,oe.M

    t0 = t0*86400.0 #convert to seconds
    t = t*86400.0

    n = sqrt(mu/a^3)
    p = a*(1-e^2)

    Wdot = (-3/2)*J2*(req/p)^2*n*cos(i)
    wdot = (3/4)*J2*(req/p)^2*n*(5*cos(i)^2-1)

    W = Wdot*(t-t0)+W0
    w = wdot*(t-t0)+W0
    M = n*(t-t0)+M0

    #Using Newton's method to solve Keplers Equation for E Eccentric Anomaly   
    f(x) = x-e*sin(x)-M
    fp(x) = 1-e*cos(x)
    E = 0.1 #first guess
    epsilon = 1f-6
    F = 1 # just to get started
    ctr = 0    
    while abs(F) > epsilon
        ctr = ctr+1
        if ctr >= 1000            
            break            
        end
        F = f(E)
        E = E - f(E)/fp(E)
    end    

    #True Anomaly
    f = 2*atan(tan(E/2)/sqrt((1-e)/(1+e)))
    
    #Flight Path Angle
    gamma = atan(e*sin(f)/(1+e*cos(f)))

    #Norm of the radius vector
    r = a*(1-e^2)/(1+e*cos(f))

    #Norm of the velocity vector
    v = sqrt(2*mu/r-mu/a)

    x = r*(cos(f+w)*cos(W)-sin(f+w)*cos(i)*sin(W))
    y = r*(cos(f+w)*sin(W)+sin(f+w)*cos(i)*cos(W))
    z = r*(sin(f+w)*sin(i))

    vx = v*(-sin(f+w-gamma)*cos(W)-cos(f+w-gamma)*cos(i)*sin(W))
    vy = v*(-sin(f+w-gamma)*sin(W)+cos(f+w-gamma)*cos(i)*cos(W))
    vz = v*(cos(f+w-gamma)*sin(i))

    rv = DataFrame(r = [x,y,z],
                    v = [vx,vy,vz])

    return rv
end

function rv2oe(rv)
    r = rv[1:3]
    v = rv[4:6]
    h = cross(r,v)
    rmag = sqrt(r'*r)
    vmag = sqrt(v'*v)
    e = cross(v,h)/mu - r/rmag
    e = sqrt(e'*e)
    a = 1/((2/rmag)-(vmag^2/mu))
    return a,e
end


function sc(u0,tspan)
    function f(u,p,t)

        x,y,z = u[1:3]
        dx,dy,dz = u[4:6]

        r = sqrt(x^2+y^2+z^2) 
        ddx = -(mu*x/r^3)*(1+3/2*J2*(req/r)^2*(1-5*(z^2/r^2)))
        ddy = -(mu*y/r^3)*(1+3/2*J2*(req/r)^2*(1-5*(z^2/r^2)))
        ddz = -(mu*z/r^3)*(1+3/2*J2*(req/r)^2*(3-5*(z^2/r^2)))
        
        du = [dx;dy;dz;ddx;ddy;ddz]
        return du
    end
    prob = ODEProblem(f,u0,tspan)
    #tp = LinRange(tspan[1],tspan[2],500) #for animating
    sol = solve(prob, reltol = 1e-9, abstol = 1e-9)
    return sol
end
