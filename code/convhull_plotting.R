library(rgl)
library(geometry)
library(purrr)

um <- matrix(c(0.68091923, 0.73223907, 0.01321388, 0.00000000,
               -0.12230255,  0.09590402,  0.98784864,  0.00000000,
               0.7220740, -0.6742611,  0.1548575,  0.0000000,
               0, 0, 0, 1), ncol = 4, nrow = 4,
             byrow = TRUE)

# Function to plot convex hull
hullPlot <- 
  function(df, vars, centroids, vertices, group, pal,
           plot_points = TRUE, plot_centroid = TRUE, plot_hull = TRUE){
    
    # Open viewing device
    open3d(windowRect = c(20, 30, 800, 800), userMatrix = um)
    palette(pal)
    
    # Plot points
    if (plot_points) {
      plot3d(x = df[,vars[1]], 
             y = df[,vars[2]], 
             z = df[,vars[3]],
             type = "s", 
             box = TRUE, 
             radius = 2, 
             col = 1, 
             xlab = "", 
             ylab = "", 
             zlab = "")
    }
    
    if (plot_centroid) {
      points3d(x = centroids[,"x"],
               y = centroids[,"y"],
               z = centroids[,"z"], 
               size = 15, 
               alpha = 0.8,
               col = centroids[,group],
               add = TRUE)
    }
    
    if (plot_hull) {
      rgl.triangles(x = vertices[,"x"],
                    y = vertices[,"y"],
                    z = vertices[,"z"],
                    col = vertices[,group],
                    alpha = 0.3, 
                    add = TRUE)
    }
    axes3d(col = "black", 
           alpha = 1)
    # # Slow down the plotting
    # Sys.sleep(0.5)
    # # Add legend
    # legend3d("top",
    #          horiz = TRUE,
    #          legend = gsub("_", " ", groups),
    #          col = pal,
    #          pch = 16,
    #          inset = c(0.02),
    #          bty = "n")
  }
