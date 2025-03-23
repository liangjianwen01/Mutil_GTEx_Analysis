``

### function to Perform QUASI ROBUST Analysis: T1 and T2
RobQuasi <- function (y,group,max2=FALSE,single.only=FALSE,...)
{
  
  y = unlist(y)
  if (max2) run = rtest2grp(y,group,...)
  else run = rtest(y,group,single.only,...)
  
  if (run$rcode==1)  return(c(run$out,run$b))
  else {
    ng = length(table(group))
    nn = ifelse(single.only,2,3)
    return(c(rep(NA,ifelse(max2,ng,nn)),run$b) ) 
  }
}
