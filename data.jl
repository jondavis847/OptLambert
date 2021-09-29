using CSV,DataFrames

Data = CSV.read(raw"C:\Users\jonda\Documents\julia\eas6722\Project\gtoc9-data\debris.csv",DataFrame)
rename!(Data,["id","t","a","e","i","W","w","M"])

function updatedebris(t)    
    D2 = copy(Data)
    
    t = t*86400.0    
    for k = 1:123        
        oe = Data[k,:]
        t0,a,e,i,W0,w0,M0 = oe.t,oe.a,oe.e,oe.i,oe.W,oe.w,oe.M
        t0 = t0*86400.0 #convert to seconds        

        mu = 398600.4418f9
        J2 = 1.08262668f-3         
        req = 6378137

        n = sqrt(mu/a^3)
        p = a*(1-e^2)        

        Wdot = (-3/2)*J2*(req/p)^2*n*cos(i)
        wdot = (3/4)*J2*(req/p)^2*n*(5*cos(i)^2-1)
    
        W = mod(Wdot*(t-t0)+W0,2*pi)
        w = mod(wdot*(t-t0)+W0,2*pi)
        M = mod(n*(t-t0)+M0,2*pi)
        D2[k,:W] = W
        D2[k,:w] = w
        D2[k,:M] = M
    end
    return D2
end

function getbyraan(low,up,t)
    D2 = updatedebris(t)
    i = (D2.W .> low) .& (D2.W .< up)
    out = D2[i,:]
    return out
end

function getwindowdata(t,D1,D2)
    tfrange = LinRange(t,t+1,72)    
    e = zeros(length(tfrange))
    m = zeros(length(tfrange))
    rp = zeros(length(tfrange))

    for i = 2:length(tfrange)
        print("...Working $i...")
        v,e[i],m[i],prop = dolambert(D1,D2,t,tfrange[i])
        a,ec = rv2oe(prop.u[end])
        rp[i] = a*(1-ec)
        print("Done\n")
    end

    Data2 = updatedebris(t)
    dW = Data2[D2+1,:W] - Data2[D1+1,:W]
    dw = Data2[D2+1,:w] - Data2[D1+1,:w]

    em = DataFrame(tf = tfrange[2:end],
                    e = e[2:end],
                    m = m[2:end],
                    rp = rp[2:end])
    return em,dW
end    

function comparedebris(D,t)
    eM = Vector{DataFrame}(undef,123)
    dW = Vector{Float64}(undef, 123)
    dw = Vector{Float64}(undef, 123)

    Threads.@threads for  i = 1:123
        print("Comparing $D to $(i-1)...\n")
        if D != i-1            
            eM[i], dW[i] ,dw[i]= getwindowdata(t,D,i-1)            
        end
    end
    return eM,dW,dw
end