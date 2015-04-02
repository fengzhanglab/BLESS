mle_estimation <- function(R,q,n)
{
  #input R is total number of reads analyzed
  #q is fraction of reads in negative control that are counted as having at
  #least one indel
  #n is number of reads in data-set that have at least one indel
  #returns MLE-estimate of cutting-frequency (between 0 and 1)

  if (R<n) 
  {
    print('Error, #indel-reads exceeds total number of reads.'); 
  }

  pspace <- (0:n)/R;
  likelihood <- array(0,n+1);

  #loop through pspace
  for (i in 1:(n+1))
  {
    #function round() is used in case round-off error above has converted
    #integers to floating-point
    likelihood[i] <- dbinom(round(n-(R*pspace[i])),round(R*(1-pspace[i])),q);    
  }
  if (sum(likelihood)>0)
  {
    likelihood <- likelihood/sum(likelihood);
    max_prob <- max(likelihood);
    mle_index <- which.max(likelihood);

    mle <- pspace[mle_index];
    std <- sqrt(sum(likelihood*(pspace^2)) - (sum(likelihood*pspace))^2);

    return(c(mle,std))    
  } else
  {
    return(c(0,0))
  }

}

