#packages required 
# uncomment to install the packages
# install.packages("gtools")
# install.packages("ggplot2")
# install.packages("gganimate")
# install.packages("magick")
#install.packages("gridExtra")




library(gtools)
library(ggplot2)
library(gganimate)
library(magick)
library(gridExtra)



plot_Gr <- function(){}
sendOutside <-function(ploter){
  plot_Gr <<- ploter
  
  
}
calcDist <- function(fVec){
  
  stV <- c(0,0)
  dist <- stV - fVec
  diffr <- sqrt(dist[1]^2 + dist[2]^2)
  
  return(round(diffr,3))
  
}

generate_transitionMat <- function(){
  
  sam_mat <- matrix(ncol = 2)
  samV <- sample(1:4,100, replace = T)
  j <- 1
  for(i in 1:100){
    j <- j+1
    sam_mat <- rbind(sam_mat,c(samV[i],samV[j])) 
    
  }
  perm_matrix <- permutations(4,2,repeats.allowed = TRUE)
  
  
  v1 <- (c(1:16))
  k_counter <- 1
  for(i in 1:16){
    prob_counter <- 0
    
    
    for(x in 1:100){
      if(identical(perm_matrix[i,],sam_mat[x,])){prob_counter <- prob_counter+1}
      
      
    }
    v1[k_counter] <- prob_counter
    k_counter <- k_counter+1
    prob_counter <- 0
  }
  
  
  pr_transition <- matrix(v1, nrow = 4, ncol = 4, byrow = TRUE ,dimnames = list(c("N","E","S","W"),c("N","E","S","W")))
  
  
  calcProb <- function(vec,sum){
    for( i in 1:4){
      vec[i] <- round(vec[i]/sum,3) 
    }
    return(vec)
  }
  
  test_mat <- matrix(c(0:0),nrow =4 , ncol = 4);
  sumVec <- c(1:4)
  for(i in 1:4){
    
    sum_row <- sum(pr_transition[i,])
    Probvec_ret <- calcProb(pr_transition[i,],sum(pr_transition[i,]))
    test_mat[i,] <- Probvec_ret
    sumVec[i] <- sum_row
    sum_row<-0
    
  }
  return(test_mat)
}


generateSteps <- function(steps = 100){
  N <- steps
  rand_var_Mat <- matrix(c(0:0),ncol = 1, nrow = N)
  i_w_mat <- matrix(c(0.2,0.4,0.3,0.1),ncol  = 4)
  trans_mat <- generate_transitionMat()
  
  for(i in 1:N){
    
    rand_var_Mat[i] <- match(1,rmultinom(1,1,i_w_mat))
    i_w_mat <- trans_mat[rand_var_Mat[i],]
    
  }
  
  return (c(rand_var_Mat))
}




calculate_Steps <- function(n.Times = 1,steps = 10){
  disp_matrix <- matrix(c(0:0),nrow = n.Times,ncol = 3)
  
  dest_path <- matrix(ncol = 4)
  for(x in 1:n.Times){
    asteps <- steps + 1
    rand.w <- generateSteps(asteps)
    
    counter <- 2
    path_plot <- matrix(c(0:0), nrow = steps+2, ncol = 2, byrow = TRUE)
    dis_v <- matrix(c(0:0), nrow = steps+2, ncol = 1, byrow = TRUE)
    steps_v <- matrix(c(0:0), nrow = steps + 2, ncol = 1,byrow = TRUE)
    
    for(i in 1:asteps){
      
      if(rand.w[i] == 2){ path_plot[i+1,] <- path_plot[i,] + c(1,0) }
      if(rand.w[i] == 3){ path_plot[i+1,] <- path_plot[i,] - c(0,1)}
      if(rand.w[i] == 1){ path_plot[i+1,] <- path_plot[i,] + c(0,1)}
      if(rand.w[i] == 4){ path_plot[i+1,] <- path_plot[i,] - c(1,0)}
      c_vec <- c(path_plot[i,])
      dist <- calcDist(c_vec)
      dis_v[i] <- dist 
      steps_v[counter] <- i
      counter <- counter + 1
      
      
    }
    disp_matrix[x,1] <- x 
    disp_matrix[x,2] <- round((sum(dis_v)/steps),3)
    disp_matrix[x,3] <- max(dis_v)
    
    path_plot <- cbind(steps_v,path_plot,dis_v)
    
    dest_path <- rbind(dest_path,path_plot)
    
    
    
    
  }
  colnames(disp_matrix) <- c("walk_number","average_displacement","largest_displacement")
  disp_matrix <- as.data.frame(disp_matrix)
  
  
  colnames(dest_path) <- c("nSteps","x_coordinate","y_coordinate","displacement")
  id <- rep(1:n.Times, each = steps + 2)
  
  dest_path <- as.data.frame(dest_path)
  
  dest_path <- cbind(dest_path[2:nrow(dest_path), ], Times =  factor(id))
  
  dest_path <- na.omit(dest_path)
  p <- ggplot(dest_path, aes (x_coordinate, y_coordinate, color = Times))  + geom_path() + geom_point(size = 2) + transition_reveal(nSteps)
  p_gif <- animate(p, width = 500, height = 500)
  
  q <- ggplot(dest_path,aes(nSteps,displacement , color = Times)) + geom_line() + geom_point(size = 4) + transition_reveal(nSteps)
  
  q_gif <- animate(q, width = 500 , height = 500)
  
  p_mgif <- image_read(p_gif)
  q_mgif <- image_read(q_gif)
  new_gif <- image_append(c(p_mgif[1], q_mgif[1]))
  for(i in 2:100){
    combined <- image_append(c(p_mgif[i], q_mgif[i]))
    new_gif <- c(new_gif, combined)
  }
  
  print(new_gif)
  
  
  sendOutside (function(){
    
    
    p <- ggplot(dest_path, aes(x = x_coordinate, y = y_coordinate, color = Times)) 
    p <- p + geom_path() + labs(x="x-coordinate in unit length", y = "y-cordinate in unit length" , title = "the position of random walkers simulated (times) times")
    
    q <- ggplot(dest_path,aes(nSteps,displacement , color = Times)) + geom_line() + labs(x = "displacement at nStep ", y = "number of steps", title = "the displacement from origin at each  step")
    
    grid.arrange(p,q,nrow = 1)
    av <- ggplot(disp_matrix,aes (x = walk_number , y = average_displacement , color = walk_number)) + geom_line(color="#69b3a2") + geom_point(size = 3 , color = "green") 
    av <-  av + labs( title = "average displacement in each walk")
    
    
    hd <- ggplot(disp_matrix,aes (x = walk_number , y = largest_displacement, color = walk_number)) + geom_line(color="#69b3a2") + geom_point(size = 3, color = "red")  
    
    hd <-  hd + labs(title = "highest displacement in each walk")
    
    grid.arrange(p,q,av,hd )
    
    
  })
  
}

calculate_Steps(10,1000)

#call the plot_Gr() function in console during simulation to render the graph 
# plot_Gr()



