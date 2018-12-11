# modify the cytometree function

CytomeTree <- function(M, minleaf = 1, t = .1, first_markers_order, verbose = TRUE)
{
  start_time = proc.time()
  
  
  #first_markers_order = vector of the the order of markers in
  # the origine matrix (from FCS file)
  
  
  if((class(M) != "matrix") & (class(M) != "data.frame"))
  {
    stop("M should be of class matrix or data.frame.")
  }
  n <- nrow(M)
  if(minleaf >= n)
  {
    stop("minleaf is superior to n.")
  }
  p <- ncol(M)
  if(p > n){
    stop("p is superior to n.")
  }
  if(any(is.na(M)))
  {
    stop("M contains NAs.")
  }
  
  if (is.null(first_markers_order))
    
  { 
  BT <- BinaryTree(M, floor(minleaf), t, verbose)
  
  annotation <- TreeAnnot(BT$labels, BT$combinations)
  
  Tree <- list("M" = M, "labels" = BT$labels,
               "pl_list"= BT$pl_list, "t"= t,
               "mark_tree" = BT$mark_tree,
               "annotation" = annotation)
  class(Tree) <- "CytomeTree"
  Tree
  }
  
  else
  {
    
  #  new code
    
  #first_markers_order = c(2,3)

  p0 = length(first_markers_order)
  p = ncol(M)
  n = nrow(M)

   # update M

  p1 = c (first_markers_order , (1:p)[!(1:p) %in% first_markers_order] ) 
  M = M[,p1]  #   new matrix  M with order : fisrt markers, others
  
   
   
   col_names <- colnames(M)
   colnames(combinations) <- col_names
   
   
   root <- tree <- mark_tree <- marks_left <- rootmarks <- list()
   
   labels <- rep(0, n)
   
   label_counter <- label_graph <- level <- 1
   # level 1 
   
   root[[level]] = [M,j]
   tree[[level]] = root
   
   
  # firt marker  = root -level 1
  # second marker = level 2
  #...
   # then algorithm BinaryTree with CytEM 
  
   
  for j in (1:p0)
  {
  
  M_j= M[,j]
  
  mc_uni <- Mclust(M_j, G=1, verbose = FALSE)
  mc_mix <- Mclust(M_j, G=2, modelNames = "E", verbose = FALSE)
  ind1 <- mc_mix$classification == 1
  ind2 <- mc_mix$classification == 2
  

  
  
  M1 <- M_j[ind1]
  M2 <- M_j[ind2]
  aic_uni <- 2*mc_uni$df - 2*mc_uni$loglik
  aic_mix <- 2*mc_mix$df - 2*mc_mix$loglik
  aic_norm_new <- (aic_uni - aic_mix)/n
  
  label <- mc_mix$classification
  mean_M1 <- mean(M1)
  mean_M2 <- mean(M2)
  pi_M1 <- length(M1)/n
  pi_M2 <- 1 - pi_M1
  if(mean_M1 > mean_M2)
    
  {
    label[ind1] <- 1
    label[ind2] <- 0
    temparameters <- c(aic_norm_new, j, mean_M2, mean_M1,
                       stats::var(M2), stats::var(M1), pi_M2, pi_M1)
  }   else   {
    label[ind1] <- 0
    label[ind2] <- 1
    temparameters <- c(aic_norm_new, j, mean_M1, mean_M2,
                       stats::var(M1), stats::var(M2), pi_M1, pi_M2)
  }
  
  child$L <-  indices[label == 0]
  child$R <-  indices[label == 1]
  
  parameters <- rbind(parameters, temparameters)

  
  }
}

nnrowpara <- nrow(parameters)
    
if(is.null(nnrowpara))
{
  return(list("mark_not_dis"= mark_not_dis))
}
if(nnrowpara > 1)
{
  parameters <- parameters[order(parameters[,1],decreasing = TRUE),]
}
return(list("mark_not_dis" = mark_not_dis, "child" = child,
            "nAIC" = parameters[,1], "ind" = parameters[,2],
            "mu1"= parameters[1,3], "mu2"= parameters[1,4],
            "Var1" = parameters[1,5], "Var2" = parameters[1,6],
            "pi1"= parameters[1,7], "pi2" = parameters[1,8]))
}














 proc.time() - start_time
}




