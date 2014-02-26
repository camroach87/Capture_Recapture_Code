# Original code courtesy of Richard Huggins, University of Melbourne
# Constructs population matrices for open populations

Pop.Mat.0 <- function(num,t,beta)
{
  
  entry <- sample(1:t,num,replace=TRUE); #uniform entry
  surv <- rgeom(num,beta);#geometric survival
  leave <- entry+surv;
  leave <- pmin(leave,t);
  
  Pop <- matrix(0,num,t);
  
  for(k in 1:num)
  {
    Pop[k,entry[k]:leave[k]] <- 1;
  }
  return(Pop);
}




Pop.Mat <- function(num,t,R,beta)
{
  
  Pop.0 <- Pop.Mat.0(num,t+R,beta);
  
  Pop <- Pop.0[,-(1:R)];#remove initial generations so stable
  
  pp <- apply(Pop,1,sum);
  
  Pop <- Pop[pp>0,];
  
  return(Pop);
}