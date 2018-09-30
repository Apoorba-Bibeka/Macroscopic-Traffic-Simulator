---
title: "Macroscopic Simulator"
author: "Apoorba Bibeka"
date: "April 1, 2017"
output:
  word_document: default
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# 0 Housekeeping 
```{r}
 rm(list=ls())
 ls()
```


# 1 Creating a macroscopic traffic simulator (LWR Model) n lane

Following are assumptions made:  
              1.Reaction time = 2 seconds  
              2.Free flow speed = 60 mph (26.8224 m/s)  
              3.$\Delta t$ = 2.5 seconds or 0.0006944444 hr
              4.$\Delta x$ = 80.4672 meters  ($\geq \Delta t \times u_f$ = 67.057~80.4672 meters or 0.05 miles)   
              
  
## 1.1 Basic macroscopic simulator functions
```{r}
macrosimulator_fun<-function(k_bumper,u_f,k_jam,del_x,del_t,len,q_input,no_iter,n
                                                          ,lane_drop=c(0,0)){
  
  #*************************************************************************************
  # # #Comment out the following values when runing the fun
  # u_f=60
  # k_jam=120
  # k_bumper=264   #Assuming average vehicle length of 20 feets
  # no_iter=5000  # Number of times steps
  # del_x=0.05
  # del_t =0.000694444
  # len= 10
  # q_input=1500
  # n=3
  #  lane_drop= c(5,2)  #At 5 miles drop 2 lane
  #*************************************************************************************
  q_initial=0
  u_initial=0
  k_initial=0
  no_sections= floor(len/del_x)
  if(lane_drop[2]!=0){
          pos_lane_drop<-lane_drop[1]/del_x
          #cell 100
          lane_drop_cell<-ceiling(pos_lane_drop)
  }
  
  Q=matrix(data=NA,ncol=no_sections,nrow=no_iter)
  q=matrix(data=NA,ncol=no_sections,nrow=no_iter)
  k=matrix(data=NA,ncol=no_sections,nrow=no_iter)
  u=matrix(data=NA,ncol=no_sections,nrow=no_iter)
  q[1,]<-q_initial
  k[1,]<-k_initial
  u[1,]<-u_initial
  s<-rep(0,length(q[1,]))
  d<-rep(0,length(q[1,]))
  
  #(no_iter-1)
  for(t in 1:(no_iter-1)){
      s= (k_bumper*n-k[t,])*del_x
    s_I_1<-2000000 # A very large no denoting unlimited capacity 
    s_i_1<-c(s[-1],s_I_1) 
    length(s_i_1)         #s_i_1 Contains the s values 
                          # for i+1 celll. Doesn't contain the value
                          #for the last cell (last cell+1)
    d= mapply(min,(q[t,]*del_t),(k[t,]*del_x))   #Find d for all cells
    #Change d for modified cell before and after lane drop
     if(lane_drop[2]!=0){
            s_i_1[(lane_drop_cell-1):(length(s_i_1)-1)]<-(k_bumper*(n-lane_drop[2])
                                      -k[t,lane_drop_cell:length(s_i_1)])*del_x
            #Last element of s_i_1 is 200,000 so do not consider only till length(s_i_1)-1
       }
       
    Q[t,]<-mapply(min, s_i_1, d)   #Get Q from d and s_i+1
    
    Q_zeroth<-q_input*del_t*n   # At each time step intialize the
                              #flow rate in the zeroth cell as q_initial
    Q_diff_1st_cell<-Q_zeroth-Q[t,1]
    Q_diff<-(-diff(Q[t,]))       #Get the flow difference of subsequent cells 

    Q_diff<-c(Q_diff_1st_cell,Q_diff)  #Get Q_diff for all the cells
    
    k[t+1,]<- k[t,] + (1/del_x)* Q_diff
      if(lane_drop[2]!=0){
    k[t+1,1:lane_drop_cell-1]<-mapply(min,k[t+1,1:(lane_drop_cell-1)],120*n)
    k[t+1,lane_drop_cell:length(k[t+1,])]<-mapply(min,k[t+1,lane_drop_cell:length(k[t+1,])],120)
    
  }
    u[t+1,]<-u_f*(1-(k[t+1,]/(n*k_jam)))
    u[t+1,]<-mapply(max,u[t+1,],0)
    if(lane_drop[2]!=0){
      u[t+1,(lane_drop_cell):length(u[t+1,])]<-u_f*(1-(k[t+1,(lane_drop_cell):length(u[t+1,])]/((n-lane_drop[2])*k_jam)))
    }
    q[t+1,]<-u[t+1,]*k[t+1,]
    
  }
  
  write.csv(q,file="apoorb.csv")
  output<-list(k,u,q)
  output
}
```

## 1.2 Runing the macroscopic simulator
```{r}
  free_flow_speed=60
  density_jam=120
  delta_x=0.05
  delta_t=0.000694444
  length_section=10
  q_input_=1500 #vehicles/hour/ln
  iterations=5000
  no_lanes=1
  results<-macrosimulator_fun(k_bumper=264,u_f=free_flow_speed,k_jam=density_jam,del_x=delta_x,
                     del_t=delta_t,len=length_section,q_input=q_input_,
                     no_iter=iterations,n=no_lanes)
  k_sim_1<-results[[1]]
  u_sim_1<-results[[2]]
  q_sim_1<-results[[3]]
  
  no_lanes=3
  results<-macrosimulator_fun(k_bumper=264,u_f=free_flow_speed,k_jam=density_jam,del_x=delta_x,
                     del_t=delta_t,len=length_section,q_input=q_input_,
                     no_iter=iterations,n=no_lanes)
  #3 lanes   
  k_sim_2<-results[[1]]
  u_sim_2<-results[[2]]
  q_sim_2<-results[[3]]
  
  
  #3 lanes with 2 lane drop at 5 miles
  no_lanes=3
  lane_drop=c(5,2)
    results<-macrosimulator_fun(k_bumper=264,u_f=free_flow_speed,k_jam=density_jam,del_x=delta_x,
                     del_t=delta_t,len=length_section,q_input=q_input_,
                     no_iter=iterations,n=no_lanes,lane_drop = lane_drop)
      k_sim_3<-results[[1]]
      u_sim_3<-results[[2]]
      q_sim_3<-results[[3]]
```

## 1.3 Function for plotting fundamental diagrams 
```{r}
plot_fundamental<-function(cell_no,sec,k_sim,q_sim,u_sim,n_lanes){
  par(mfrow=c(1,1))
  #Flow Density Curve  
  plot(k_sim[-1,cell_no],q_sim[-1,cell_no],xlab="Density(Veh/Mi/ln)",
         ylab="Flow Rate (Veh/Hr/ln)"
         ,main="Flow vs. Density Curve",
          sub=paste("section =",sec,"miles",n_lanes,sep=" " ),cex.sub=0.75,
        type="p",col="dark blue")
         lines(k_sim[-1,cell_no],q_sim[-1,cell_no],col="dark red",lwd="2")
  #Speed Density Curve
  plot(k_sim[-1,cell_no],u_sim[-1,cell_no],xlab="Density(Veh/Mi/ln)",
         ylab="Speed (Mph)"
         ,main="Speed vs. Density Curve",
          sub=paste("section =",sec,"miles",n_lanes,sep=" " ),cex.sub=0.75,
          type="p",col="dark blue")
        lines(k_sim[-1,cell_no],u_sim[-1,cell_no],col="dark red",lwd="2")

  #Speed Flow Curve
  plot(q_sim[-1,cell_no],u_sim[-1,cell_no],xlab="Flow Rate (Veh/Hr/ln)",
         ylab="Speed (Mph)"
         ,main="Speed vs. Flow Rate Curve",
        sub=paste("section =",sec,"miles",n_lanes,sep=" " ),cex.sub=0.75,
         type="p",col="dark blue")
        lines(q_sim[-1,cell_no],u_sim[-1,cell_no],col="dark red",lwd="2")
}
```


## 1.4 Fundamental Diagram at x=4.5, 5.0 and 5.5 miles for 1 lane roadway   

From the following figure it can be seen that the fundamental diagrams stays the same for 4.5,  
5.0 and 5.5 miles section.
```{r}
  #Find cell numbers which 4.5, 5 and 5.5 miles section
  cell_4_5<-ceiling(4.5/delta_x)
  cell_5<-ceiling(5/delta_x)
  cell_5_5<-ceiling(5.5/delta_x)
  
  #Ploting the fundamental diagrams for 1 lanes
  plot_fundamental(cell_4_5,4.5,k_sim_1,q_sim_1,u_sim_1,"1 lanes")
  plot_fundamental(cell_5,5,k_sim_1,q_sim_1,u_sim_1,"1 lanes")
  plot_fundamental(cell_5_5,5.5,k_sim_1,q_sim_1,u_sim_1,"1 lanes")
```

# 2 Fundamental Diagram at x=4.5, 5.0 and 5.5 miles for 3 lane roadway   

From the following figure it can be seen that the fundamental diagram stays the same for 4.5,  
5.0 and 5.5 miles section. The fundamental digrams for the 3 lane case and 1 lane case are also the same.
```{r}
    #Ploting the fundamental diagrams for 3 lanes
  plot_fundamental(cell_4_5,4.5,k_sim_2/3,q_sim_2/3,u_sim_2,"3 lanes")
  plot_fundamental(cell_5,5,k_sim_2/3,q_sim_2/3,u_sim_2,"3 lanes")
  plot_fundamental(cell_5_5,5.5,k_sim_2/3,q_sim_2/3,u_sim_2,"3 lanes")
  
  k<-k_sim_2/3
```


# 3 Fundamental Diagram at x=4.5, 5.0 and 5.5 miles for lane drop case (3 lane roadway with 2 lane drop at 5 miles)   

For 4.5 miles section the fundamental diagrams is same as the fundamental diagram obtained for the 3 lane case.  
For 5 miles section the flow increases to the maximum flow in q-k and u-q curve and then it drops to zero.  
For 5.5 miles section the fundamental diagrams is like that for the 1 lane scenario. The only difference is that the maximum flow rate for this case is more than 1500 veh/hr/ln.
```{r}
     pos_lane_drop<-lane_drop[1]/delta_x
          #cell 100
          lane_drop_cell<-ceiling(pos_lane_drop)
  k_sim_3[,1:(lane_drop_cell-1)]=k_sim_3[,1:(lane_drop_cell-1)]/3
  q_sim_3[,1:(lane_drop_cell-1)]=q_sim_3[,1:(lane_drop_cell-1)]/3
  
  plot_fundamental(cell_4_5,4.5,k_sim_3,q_sim_3,u_sim_3,"3 lanes (no drop)")
  plot_fundamental(cell_5,5,k_sim_3,aa,u_sim_3,"1 lane (Lane drop cell)")
  plot_fundamental(cell_5_5,5.5,k_sim_3,q_sim_3,u_sim_3,"1 Lane section")
```


