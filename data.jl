using CSV,DataFrames

Data = CSV.read(raw"C:\Users\jonda\Documents\julia\eas6722\Project\gtoc9-data\debris.csv",DataFrame)
rename!(Data,["id","t","a","e","i","W","w","M"])
insertcols!(Data,9, :state => trues(123))
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

function getortho(t,D)
    D2 = updatedebris(t)
    e = 0.25
    poslow = mod(D2[D+1,:W] + pi/2 - e, 2*pi)    
    posup = mod(D2[D+1,:W] + pi/2 + e, 2*pi)    
    neglow = mod(D2[D+1,:W] - pi/2 - e, 2*pi)
    negup = mod(D2[D+1,:W] - pi/2 + e, 2*pi)
    
    i = ((D2.W .> poslow) .& (D2.W .< posup) .| (D2.W .> neglow) .& (D2.W .< negup)) .& D2.state
    # if we dont find any, find the most ortho debris
    if sum(i) < 1        
        W = D2[D+1,:W]
        delW = abs.(D2.W.-W)
        mini = min(delW[D2.state]...)
        i = delW .== mini
    end    
    out = D2[i,:]
    return out    
end

function getwindowdata(t,D1,D2,n)
    tfrange = LinRange(t+0.5,t+0.7,n)    
    e = zeros(length(tfrange))
    m = zeros(length(tfrange))
    rp = zeros(length(tfrange))

    for i = 1:length(tfrange)        
        v,e[i],m[i],prop = dolambert(D1,D2,t,tfrange[i])
        a,ec = rv2oe(prop.u[end])
        rp[i] = a*(1-ec)        
    end

    #Data2 = updatedebris(t)
    #dW = Data2[D2+1,:W] - Data2[D1+1,:W]
    #dw = Data2[D2+1,:w] - Data2[D1+1,:w]
    id = D2*ones(length(tfrange))
    em = DataFrame(tf = tfrange,
                    e = e,
                    m = m,
                    rp = rp,
                    id = id)

    return em#,dW,dw
end    

function comparedebris(D,t)    
    O = getortho(t,D)          
    n = size(O,1)
    if n < 4
        m = 25
    else
        m = 5
    end
    #print("Found $n debris with roughly orthogonal orbits\n")
    eM = Vector{DataFrame}(undef,size(O,1))    

    Threads.@threads for  i = 1:size(O,1)
        #print("Comparing $D to $(i-1)...\n")        
        eM[i] = getwindowdata(t,D,O[i,:].id,m)
    end
    return eM
end

function getraanstats(dW,EM)
    m = zeros(123)
    sd = zeros(123)
    mini = zeros(123)
    maxi = zeros(123)
    s3p = zeros(123)
    s3m = zeros(123)
    
    for i = 2:123
        validerror = EM[i].e .< 100
        validrp = EM[i].rp .> 6.6e6
        v = (validerror .& validrp)
        m[i] = mean(EM[i].m[v])
        sd[i] = std(EM[i].m[v])
        mini[i] = min(EM[i].m[v]...)
        maxi[i] = max(EM[i].m[v]...)        
    end   

    scatter(dW[2:end],m[2:end], xlabel = "RAAN diff (rad)", ylabel = "(kg)", title = "Piece 0 to all, T = 24000 days, 1 day window, 20 min increment", label = "mean", size = (1000,600))    
    scatter!(dW[2:end],mini[2:end],color = :orange, label = "min")
    scatter!(dW[2:end],maxi[2:end],color = :red, label = "max")
    
end

function gettimestats(t0,EM)
    n = size(EM,1)
    m = zeros(n)    
    mini = zeros(n)
    maxi = zeros(n)
    tf = LinRange(t0+0.5,t0+.7,5) 
    tf = tf[2:end]
    for i = 1:n
        validerror = EM[i].e .< 100
        validrp = EM[i].rp .> 6.6e6
        v = (validerror .& validrp)
        m[i] = mean(EM[i].m[v])        
        mini[i] = min(EM[i].m[v]...)
        maxi[i] = max(EM[i].m[v]...)        
    end   
    scatter(tf[2:end],m[2:end], xlabel = "Tf (days)", ylabel = "(kg)", title = "Piece 0 to ortho, T = 24000 days, 0.2 day window, 60 min increment", label = "mean", size = (1000,600))    
    scatter!(tf[2:end],mini[2:end],color = :orange, label = "min")
    scatter!(tf[2:end],maxi[2:end],color = :red, label = "max")
end