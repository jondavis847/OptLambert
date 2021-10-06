using ProgressBars

function getnextpiece(D,t)    
    comp = comparedebris(D,t)
    mins = map(x -> min(x.m...),comp)
    mini = mins .== min(mins...)
    i = comp[mini][1].m .== min(mins...)
    nextpiece = convert(Int64,comp[mini][1].id[i][1])
    arrivetime = comp[mini][1].tf[i][1]
    cost = comp[mini][1].m[i][1]
    return nextpiece,arrivetime,cost
end

function planmission(cp,t,N,m0)
    debrisleft = sum(Data.state)
    if (debrisleft > N) & (debrisleft < 2*N)
        N = round(debrisleft/2)
    elseif debrisleft < N
        N = debrisleft
    end

    N = convert(Int64,N)

    at = Vector{Float64}(undef,N)
    dt = Vector{Float64}(undef,N)
    id = Vector{Int64}(undef,N)
    c = Vector{Float64}(undef,N)
    m = Vector{Float64}(undef,N)
    msc = Vector{Float64}(undef,N)
    #sol = Vector{OrdinaryDiffEq.ODECompositeSolution}(undef,N)

    at[1] = t
    id[1] = cp
    msc[1] = m0
    c[1] = 0
    m[1] = 0
    dt[1] = at[1]+5
    Data[cp+1,:].state = false
    broke = false
    
    for de in ProgressBar(1:N-1)    
        ctr = de
        id[de+1],at[de+1],c[de+1] = getnextpiece(id[de],dt[de])
        if msc[de]-c[de+1] < 2000
            broke = true
            break
        end
        #s,e,m[de+1],sol[de+1] = dolambert(id[de],id[de+1],dt[de],at[de+1])
        s,e,m[de+1] = dolambert(id[de],id[de+1],dt[de],at[de+1])
        msc[de+1] = msc[de] - 30 - m[de+1]    
        Data[id[de+1]+1,:].state = false
        dt[de+1] = at[de+1]+5        
    end

    J = 55+alpha*(msc[1]-mdry)^2
    o = DataFrame(id=id,at=at,dt=dt,c=c,m=m,msc=msc)#,sol=sol)

    if broke
        delete!(o,ctr+1:size(o,1))
    end
    return o,J
end

function reset()
    Data.state = trues(123)
end

function fullrun(dn,dp)
    reset()
    N = 11
    delW = pi/3    
    t0 = 23467
    mis = Vector{DataFrame}(undef,N)
    J = Vector{Float64}(undef,N)
    
    
    for i in ProgressBar(1:N)
        n = 12-dn[i]
        m0 = mdry+n*mde+(mp-dp[i]+5)

        nextpiece = Data[Data.state,:][end,:]
        mis[i],J[i] = planmission(nextpiece.id,t0,n,m0)
        t0 = mis[i].dt[end] + 30
    end
    return mis,J
end

        