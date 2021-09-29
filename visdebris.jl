using Plots

function visdebris(sc)    
    t = sc.t

    #t = LinRange(24000.0,24001.0,500)    
    p = zeros(3,size(Data)[1], length(t))
    for i = 1:length(t)
        for d = 1:size(Data)[1]
            p[:,d,i]=oe2rv(t[i],Data[d,:]).r
        end
    end

    deb = [1+1,42+1]
    anim = @animate for i = 1:length(t)

        plot(
            background_color = "black",         
            title = string(round(t[i]/86400, digits = 3)),
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
         
        scatter!(p[1,deb,i],p[2,deb,i],
            markercolor = "yellow",                     
            markersize = 5,
            markerstrokewidth = 0,                
            subplot = 1
        )

        scatter!([sc.u[i][1]],[sc.u[i][2]],
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
         
        scatter!(p[1,deb,i],p[3,deb,i],
            markercolor = "yellow",         
            markerstrokewidth = 0,                 
            markersize = 5,            
            subplot = 2
        )

        scatter!([sc.u[i][1]],[sc.u[i][3]],
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
         
        scatter!(p[2,deb,i],p[3,deb,i],
            markercolor = "yellow",                     
            markerstrokewidth = 0,              
            markersize = 5,
            subplot = 3
        )

        scatter!([sc.u[i][2]],[sc.u[i][3]],
            markercolor = "lime",
            markersize = 5,
            markerstrokewidth = 0,                          
            subplot = 3
        )

    end

    gif(anim,"visdebris.gif",fps = 30)
    
end
