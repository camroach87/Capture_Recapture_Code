Pop.Mat.0 <- function(num,T,beta)
{
  
  entry <- sample(1:T,num,replace=TRUE); #uniform entry
  
  surv <- rgeom(num,beta);#geometric survival
  
  leave <- entry+surv;
  
  leave <- pmin(leave,T);
  
  Pop <- matrix(0,num,T);
  
  for(k in 1:num)
  {
    Pop[k,entry[k]:leave[k]] <- 1;
  }
  Pop;
}




Pop.Mat <- function(num,T,R,beta)
{
  
  Pop.0 <- Pop.Mat.0(num,T+R,beta);
  
  Pop <- Pop.0[,-(1:R)];#remove initial generations so stable
  
  pp <- apply(Pop,1,sum);
  
  Pop <- Pop[pp>0,];
  
  Pop;
}

####################################################

beta <- 0.01;

num <- 10000;

T <- 1000;

R <- 500;

Pop <- Pop.Mat(num,T,R,beta);

N <- apply(Pop,2,sum);

plot(table(apply(Pop,1,sum)))


plot(N,type="l");
