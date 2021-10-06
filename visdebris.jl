using Plots

function visdebris(mis)   
    t = Vector{Float64}(undef,0)
    deb = Matrix{Int64}(undef,(0,2))
    scr = Matrix{Float64}(undef,(0,3))
    waittime = 5 #sec
    transfertime = 30 #sec
    fps = 30 
    for m = 1:size(mis,1)
        d1 = mis[m,:id]
        waiting = LinRange(mis[m,:at],mis[m,:dt],waittime*fps)
        for w = 1:length(waiting)
            scr = [scr; oe2rv(waiting[w],d1).r']
            deb = vcat(deb, [nothing nothing])
        end        
        t = [t;waiting]        
        if m != size(mis,1)
            transfer = LinRange(mis[m,:dt],mis[m+1,:at],transfertime*fps)
            sol = mis[m+1,:sol](transfer.*86400)
            for tr = 1:length(transfer)
                scr = [scr;sol.u[tr][1:3]']
                deb = vcat(deb, [mis[m,:id],mis[m+1,:id]]')
            end            
            t = [t;transfer]            
        end
    end
    
    p = zeros(3,size(Data)[1], length(t))
    for i = 1:length(t)
        for d = 1:size(Data)[1]            
            p[:,d,i]=oe2rv(t[i],d-1).r            
        end
    end
    return scr
    anim = @animate for i in ProgressBar(1:length(t))

        plot(
            background_color = "black",         
            title = string(round(t[i], digits = 3)),
            showaxis = false,
            grid = false, 
            ticks = false,
            legend = false,                 
            widen = false ,  
            lims = (-1f7, 1f7),         
            size = (500*3,500*1.05),
            layout = (1,3)
        )

        scatter!(p[1,:,i],p[2,:,i],              
            markersize = 3,
            marker_z = p[3,:,i],         
            markercolor = :blues,
            markerstrokewidth = 0,         
            subplot = 1      
        )         
        
        if !any(isnothing.(deb[i,:]))
            scatter!(p[1,deb[i,:].+1,i],p[2,deb[i,:].+1,i],
                markercolor = "yellow",                     
                markersize = 5,
                markerstrokewidth = 0,                
                subplot = 1
            )
        end

        scatter!([scr[i,1]],[scr[i,2]],
            markercolor = "lime",            
            markersize = 5,
            markerstrokewidth = 0,                  
            subplot = 1
        )

        scatter!(p[1,:,i],p[3,:,i],                  
            markersize = 3,
            marker_z = p[2,:,i],         
            markercolor = :blues,
            markerstrokewidth = 0,         
            subplot = 2     
        )
         
        if !any(isnothing.(deb[i,:]))
            scatter!(p[1,deb[i,:].+1,i],p[3,deb[i,:].+1,i],
                markercolor = "yellow",         
                markerstrokewidth = 0,                 
                markersize = 5,            
                subplot = 2
            )
        end
        scatter!([scr[i,1]],[scr[i,3]],
            markercolor = "lime",
            markerstrokewidth = 0,                  
            markersize = 5,            
            subplot = 2
        )

        scatter!(p[2,:,i],p[3,:,i],                  
            markersize = 3,
            marker_z = p[1,:,i],         
            markercolor = :blues,
            markerstrokewidth = 0,              
            subplot = 3     
        )
         
        if !any(isnothing.(deb[i,:]))
            scatter!(p[2,deb[i,:].+1,i],p[3,deb[i,:].+1,i],
                markercolor = "yellow",                     
                markerstrokewidth = 0,              
                markersize = 5,
                subplot = 3
            )
        end
        scatter!([scr[i,2]],[scr[i,3]],
            markercolor = "lime",
            markersize = 5,
            markerstrokewidth = 0,                          
            subplot = 3
        )

    end

    gif(anim,"visdebris.gif",fps = 30)
    
end
