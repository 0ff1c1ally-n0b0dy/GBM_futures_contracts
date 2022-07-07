#############
#Libraries
#############
library(foreach)
library(doParallel)
library(ExtDist)
set.seed(123456789)

###########
#Change the following path to the path where this code is located in your laptop. 
###########
setwd("/home/mandrei/Desktop/Futures")

################
#Initialize the parameters
##################
rf<-0.00747
C<-0.01       #commission and slippage 
M<-5
G<-5
P0<-100
Q0<-5
delta_t<-1/2530
L<-10
L_range<-seq(from=0,to=L,by=delta_t)
N<-10^(4)

r_range<-c(1,2)
len_r_range<-length(r_range)

mu_range<-seq(from=0.07,to=0.19,by=0.03)
#mu_range<-seq(from=0,to=0.25,length.out=2)
len_mu_range<-length(mu_range)

sigma_range<-seq(from=0.07,to=0.13,by=0.02)
#sigma_range<-seq(from=0.1,to=0.1,length.out=1)
len_sigma_range<-length(sigma_range)

#u_range<-seq(from=110,to=300,by=10)
u_range<-seq(from=300,to=300,length.out=1)
#u_range<-P0+(1:7)*M
len_u_range<-length(u_range)

alpha_range<-seq(from=0.1,to=0.7,by=0.02)
#alpha_range<-seq(from=0.05,to=1,length.out=2)
#alpha_range<-c(0.05,0.1,0.15)
len_alpha_range<-length(alpha_range)

n_iters<-len_mu_range*len_sigma_range*len_u_range*len_alpha_range
count_iters<-1
###############
#Specify the number of cores to use and create the progress bar. Initialization goes into a function that goes into the options of the doSNOW package.
#Too many cores and the laptop freezes up. 
##############
n_cores<-120
filename_txt<-paste0("Printing_pilot_study6_new_formulas,u=300,M=",M,",G=",G,",P0=",P0,",delta_t=",round(delta_t,digits=5),",L=",round(L,digits=3),",N=",N,".txt")
filename_err<-paste0("Error_pilot_study6_new_formulas,u=300,M=",M,",G=",G,",P0=",P0,",delta_t=",round(delta_t,digits=5),",L=",round(L,digits=3),",N=",N,".txt")
cl <- makeCluster(n_cores)
registerDoParallel(cl)
err <- file(filename_err, open="wt")
sink(err, type="message")

writeLines(c(""), filename_txt)

results<-foreach(mu=mu_range,.combine = "rbind",.export="rLaplace")%:%
  foreach(sigma=sigma_range,.combine="rbind")%:%
  foreach(u=u_range,.combine="rbind")%:%
  foreach(alpha=alpha_range,.combine="rbind")%dopar%
  {
    message<-paste0(round(100*n_cores*(count_iters-1)/n_iters,digits=4),"%.\t Started mu=",mu,", sigma=",sigma,", u=",u,", alpha=",alpha," at time ",Sys.time(),"\n")
    cat(message,file = filename_txt,append = T)
    Returns_whole<-matrix(0,nrow = N,ncol = 1)
    Returns_nwhole<-matrix(0,nrow = N,ncol = 1)
    crra_whole<-matrix(0,nrow = N,ncol = len_r_range)
    crra_nwhole<-matrix(0,nrow = N,ncol = len_r_range)
    R<-(M/alpha)-G
    omega0<-((M+alpha*C)/alpha)*Q0
    count_stop_loss<-0
    count_objective<-0
    count_end<-0
    for(i in 1:N)
    {
      P<-c()
      P<-c(P,P0)
      t<-0
      flag_obj<-1
      flag_trail<-1
      
      while(t<=L)
      {
        #z<-rnorm(n=1,mean=0,sd=1)
        z<-rLaplace(n=1,mu=0,b=delta_t)
        P<-c(P,P[length(P)]*exp(mu*delta_t+sigma*z))
        len_P<-length(P)
        if (P[len_P]<=max(P)-R)
        {
          v<-t+delta_t
          Pv<-P[len_P]
          Pv_prime<-max(P)
          flag_trail<-0
          break
        }
        if (P[len_P]>=u)
        {
          v<-t+delta_t
          Pv<-P[len_P]
          Pv_prime<-max(P)
          flag_obj<-0
          break
        }
        t<-t+delta_t 
      }
      
      #1
      if (flag_obj==1 && flag_trail==1)
      {
        v<-L
        Pv<-P[len_P]
        Pv_prime<-max(P)
        count_end<-count_end+1
      }
      #2
      #if obj is triggered first
      if (flag_obj==0 && flag_trail==1)
        count_objective<-count_objective+1    
      #3
      #if trail is triggered first
      if (flag_obj==1 && flag_trail==0)
        count_stop_loss<-count_stop_loss+1
      
      Returns_nwhole[i]<-log((M/(M+C*alpha)+(alpha/(M+C*alpha))*(Pv-Pv_prime))*exp((alpha/(M+C*alpha))*(Pv_prime-P[1])))
      w_nwhole<-exp(Returns_nwhole[i])*omega0
      #w_nwhole<-(M/(M+C*alpha)+(alpha/(M+C*alpha))*(Pv-Pv_prime))*omega0*exp((alpha/(M+C*alpha))*(Pv_prime-P[1]))
      
      a<-(alpha/(M+C*alpha))*(Pv-P[1])+log(Q0-1)+1/(2*Q0-2)
      w_whole<-(((M/alpha)-Pv_prime+Pv))*(1+exp(a-0.5*exp(-a)-0.25*exp(-2*a)))
      Returns_whole[i]<-log(w_whole/omega0)
      #Returns_whole[i]<-log((((M/alpha)-Pv_prime+Pv)/omega0)*(1+exp(a-0.5*exp(-a)-0.25*exp(-2*a))))
      
      for(j in 1:len_r_range)
      {
        if(r_range[j]==1)
        {
          crra_nwhole[i,j]<-log(w_nwhole)
          crra_whole[i,j]<-log(w_whole)
        }
        else
        {
          crra_nwhole[i,j]<-(1/(1-r_range[j]))*(w_nwhole^(1-r_range[j]))
          crra_whole[i,j]<-(1/(1-r_range[j]))*(w_whole^(1-r_range[j]))
        }
      }
    }
    
    mean_return_whole<-mean(Returns_whole)
    sd_return_whole<-sd(Returns_whole)
    sharpe_whole<-(mean_return_whole-rf)/sd_return_whole
    modigliani_whole<-sharpe_whole*sd_return_whole+rf-mean_return_whole
    
    mean_return_nwhole<-mean(Returns_nwhole)
    sd_return_nwhole<-sd(Returns_nwhole)
    sharpe_nwhole<-(mean_return_nwhole-rf)/sd_return_nwhole
    modigliani_nwhole<-sharpe_nwhole*sd_return_nwhole+rf-mean_return_nwhole
    
    count_iters<-count_iters+1 
    
    cbind(rep(mu,times=len_r_range),rep(sigma,times=len_r_range),rep(u,times=len_r_range),rep(alpha,times=len_r_range),r_range,
          rep(mean_return_nwhole,times=len_r_range),rep(mean_return_whole,times=len_r_range),
          rep(sd_return_nwhole,times=len_r_range),rep(sd_return_whole,times=len_r_range),
          apply(crra_nwhole,2,mean),apply(crra_whole,2,mean),
          rep(sharpe_nwhole,times=len_r_range),rep(sharpe_whole,times=len_r_range),
          rep(modigliani_nwhole,times=len_r_range),rep(modigliani_whole,times=len_r_range),
          rep(count_objective/N,times=len_r_range),rep(count_stop_loss/N,times=len_r_range),rep(count_end/N,times=len_r_range))
    
  }

colnames(results)<-c("mu","sigma","u","alpha","r",
                     "mean return not whole","mean return whole",
                     "sd return not whole", "sd return whole",
                     "mean crra not whole","mean crra whole",
                     "sharpe ratio not whole","sharpe ratio whole",
                     "modigliani not whole","modigliani whole",
                     "prob obj","prob trail","prob neither")

(results<-as.data.frame(results))

cat(paste0("The maximum expected CRRA for whole number of contracts is ",max(results$`mean crra whole`)),"\n")
cat(paste0("The maximum expected CRRA for not whole number of contracts is ",max(results$`mean crra not whole`),"\n"))

cat(paste0("The parameters that correspond to the maximum expected CRRA for whole number of contracts are: mu=",results$mu[which.max(results$`mean crra whole`)],
             ",sigma=",results$sigma[which.max(results$`mean crra whole`)],",u=",results$u[which.max(results$`mean crra whole`)],
             ",alpha=",results$alpha[which.max(results$`mean crra whole`)],",r=",results$r[which.max(results$`mean crra whole`)],"\n"))

cat(paste0("The parameters that correspond to the maximum expected CRRA for not whole number of contracts are: mu=",results$mu[which.max(results$`mean crra not whole`)],
             ",sigma=",results$sigma[which.max(results$`mean crra not whole`)],",u=",results$u[which.max(results$`mean crra not whole`)],
             ",alpha=",results$alpha[which.max(results$`mean crra not whole`)],",r=",results$r[which.max(results$`mean crra not whole`)],"\n"))


cat(paste0("For this maximum expected CRRA for whole number of contracts, the mean return=",results$`mean return`[which.max(results$`mean crra whole`)],
             ",sd return=",results$`sd return`[which.max(results$`mean crra whole`)],",prob obj=",results$`prob obj`[which.max(results$`mean crra whole`)],
             ",prob trail=",results$`prob trail`[which.max(results$`mean crra whole`)],",prob neither=",results$`prob neither`[which.max(results$`mean crra whole`)],"\n"))      

cat(paste0("For this maximum expected CRRA for not whole number of contracts, the mean return=",results$`mean return`[which.max(results$`mean crra not whole`)],
             ",sd return=",results$`sd return`[which.max(results$`mean crra not whole`)],",prob obj=",results$`prob obj`[which.max(results$`mean crra not whole`)],
             ",prob trail=",results$`prob trail`[which.max(results$`mean crra not whole`)],",prob neither=",results$`prob neither`[which.max(results$`mean crra not whole`)],"\n"))      

filename_csv<-paste0("Results_pilot_study6_new formulas,u=300,M=",M,",G=",G,",P0=",P0,",delta_t=",round(delta_t,digits=5),",L=",round(L,digits=3),",N=",N,".csv")
write.csv(results,filename_csv)

sink(type="message")
close(err)
#stop the cluster of cores assigned for this task
stopCluster(cl) 

