## Generate simulted data

#### INT - DGP 1
#### SLOPE - DGP 2
#### INTESLOPE - DGP3 

ok=FALSE
while(!ok){
  gr=list()
  caso='inteslope'
  gr[[1]]=rep("A",100)
  gr[[2]]=rep("B",100)
  gr[[3]]=rep("C",100)
  gr[[4]]=rep("D",90)
  gr[[5]]=rep("E",90)
  gr[[6]]=rep("F",90)
  gr[[7]]=rep("G",90)
  gr[[8]]=rep("H",95)
  gr[[9]]=rep("I",95)
  gr[[10]]=rep("L",95)
  
  x = rnorm(100, mean = 0, sd = 1)
  x=sort(x)
  
  f = rnorm(100, mean = 0, sd = 1)
  w = rnorm(90, mean = 0, sd = 1)
  w=sort(w)
  
  g= rnorm(90, mean = 0, sd = 1)
  
  a = rnorm(95, mean=0, sd =1)
  
  a= sort(a)
  h = rnorm(95, mean=0, sd=1)
  
  
  if(caso=='int'){
    b = -7 + 4*x - 3*f
    r = -7  + 4*x - 3*f
    z = -7  + 4*x - 3*f
    n = -4 + 4*w - 3*g
    m = -4 +4*w - 3*g
    o = -4 + 4*w - 3*g
    q = -4 + 4*w - 3*g
    c = -2  +4*a - 3*h
    d = -2 +4*a - 3*h
    e = -2 +4*a - 3*h
    ## secondo confronto, definisco gli eta
    b2 = -5 - 2*x + 2*f
    r2 = -5 - 2*x + 2*f
    z2 = -5 - 2*x + 2*f
    n2 = -5 - 2*w + 2*g
    m2 = -5 - 2*w + 2*g
    o2 = -5 - 2*w + 2*g
    q2 = -5 - 2*w + 2*g
    c2 = -2 -2*a + 2*h
    d2 = -2 -2*a + 2*h
    e2 = -2 -2*a + 2*h
  }
  
  if(caso=='slope'){
    b = -1 + 5*x - 3*f #+ rnorm(100, 0, 0.02)
    r = -1  + 5*x - 3*f #+ rnorm(100, 0, 0.02)
    z = -1  + 5*x - 3*f #+ rnorm(100, 0, 0.02)
    n = -1 + 2*w - 3*g #+ rnorm(90, 0, 0.02)
    m = -1 +2*w - 3*g  #+ rnorm(90, 0, 0.02)
    o = -1 + 2*w - 3*g #+ rnorm(90, 0, 0.02)
    q = -1 + 2*w - 3*g #+ rnorm(90, 0, 0.02)
    c = -1  -1*a - 3*h #+ rnorm(95, 0, 0.02)
    d = -1  -1*a - 3*h #+ rnorm(95, 0, 0.02)
    e = -1  -1*a - 3*h #+ rnorm(95, 0, 0.02)
    ## secondo confronto, definisco gli eta
    b2 = -2 - 2*x + 2*f #+ rnorm(100, 0, 0.02)
    r2 = -2 - 2*x + 2*f #+ rnorm(100, 0, 0.02)
    z2 = -2 - 2*x + 2*f #+ rnorm(100, 0, 0.02)
    n2 = -2 - 2*w + 2*g #+ rnorm(90, 0, 0.02)
    m2 = -2 - 2*w + 2*g #+ rnorm(90, 0, 0.02)
    o2 = -2 - 2*w + 2*g #+ rnorm(90, 0, 0.02)
    q2 = -2 - 2*w + 2*g #+ rnorm(90, 0, 0.02)
    c2 = -2 -6*a + 2*h #+ rnorm(95, 0, 0.02)
    d2 = -2 -6*a + 2*h #+ rnorm(95, 0, 0.02)
    e2 = -2 -6*a + 2*h #+ rnorm(95, 0, 0.02)
    
  }
  
  if(caso=='inteslope'){
    b = -6 + 5*x - 5*f  #+ rnorm(100, 0, 0.02)
    r = -6  + 5*x - 5*f #+ rnorm(100, 0, 0.02)
    z = -6  + 5*x - 5*f #+ rnorm(100, 0, 0.02)
    n = -4 + 2*w - 5*g #+ rnorm(90, 0, 0.02)
    m = -4 + 2*w - 5*g #+ rnorm(90, 0, 0.02)
    o = -4 + 2*w - 5*g #+ rnorm(90, 0, 0.02)
    q = -4 + 2*w - 5*g #+ rnorm(90, 0, 0.02)
    c = -8  -1*a - 5*h #+ rnorm(95, 0, 0.02)
    d = -8  -1*a - 5*h #+ rnorm(95, 0, 0.02)
    e = -8  -1*a - 5*h #+ rnorm(95, 0, 0.02)
    ## secondo confronto, definisco gli eta
    b2 = +1 -4*x + 2*f #+ rnorm(100, 0, 0.02)
    r2 = +1 -4*x + 2*f #+ rnorm(100, 0, 0.02)
    z2 = +1 -4*x + 2*f #+ rnorm(100, 0, 0.02)
    n2 = +1 - 4*w + 2*g #+ rnorm(90, 0, 0.02)
    m2 = +1 - 4*w + 2*g #+ rnorm(90, 0, 0.02)
    o2 = +1 - 4*w + 2*g #+ rnorm(90, 0, 0.02)
    q2 = +1 - 4*w + 2*g #+ rnorm(90, 0, 0.02)
    c2 = -1 +2*a + 2*h #+ rnorm(95, 0, 0.02)
    d2 = -1 +2*a + 2*h #+ rnorm(95, 0, 0.02)
    e2 = -1 +2*a + 2*h #+ rnorm(95, 0, 0.02)
  }
  
  
  
  
  
  
  ## ora genero le probabilit?
  t_rude = list(x,x,x,w,w,w,w,a,a,a)
  t2_rude = list(f,f,f,g,g,g,g,h,h,h)
  N = length(t_rude)
  
  n_rig_gruppoj = rep(0,N)
  for(i in 1:N){
    n_rig_gruppoj[i] = length(t_rude[[i]])
  }
  
  dim_max = max(n_rig_gruppoj)
  dim_max
  
  nt_min = min(n_rig_gruppoj)
  nt_min
  
  eta1 = list(b,r,z,n,m,o,q,c,d,e)
  eta2 = list(b2,r2,z2,n2,m2,o2,q2,c2,d2,e2)
  p1=matrix(NA, nrow= 10,ncol = 100)
  p2=matrix(NA, nrow= 10,ncol = 100)
  
  for(i in 1:10){
    p1[i,1:n_rig_gruppoj[i]] = exp(eta1[[i]])/(1 + exp(eta1[[i]]) + exp(eta2[[i]])) #prob di app alla calsse 1
    p2[i,1:n_rig_gruppoj[i]] = exp(eta2[[i]])/(1 + exp(eta1[[i]]) + exp(eta2[[i]])) # // // // // // // // //2
  }
  
  p0 = 1 - (p1 + p2) #ovviamente la somma edda ess 1
  p0[p0<0]=0
  
  ## adesso genero le risposte
  K=3
  y = matrix(NA, nrow=N, ncol=dim_max)
  
  for(i in 1:N){
    for(j in 1:n_rig_gruppoj[i]){
      pout = c(p0[i,j],p1[i,j],p2[i,j])
      y[i,j] = sample(1:K, 1, pout, replace = T)
    }
  }
  col0=c()
  col1=c()
  col2=c()
  for(i in 1:10){
    col1=c(col1,t_rude[[i]])
    col2=c(col2,t2_rude[[i]])
    col0=c(col0,na.omit(y[i,]))
  }
  
  
  dataf=data.frame(col0,col1,col2,unlist(gr))
  names(dataf)<-c("y","t","x","gr")
  dataf$y<-as.character(dataf$y)
  dataf$t<-as.double(dataf$t)
  dataf$x<-as.double(dataf$x)
  for(i in 1:10){
    cat=table(y[i,])
    ok=(cat[1]>10 & cat[2]>10 & cat[3]>10)
  }
  
}

#for each case in (int, slope, inteslope) run the following code:

simulation_study<-JMSPME(dataset = dataf,resp_var_name = "y", ran_var_name = "t",groups_name = "gr",case=caso)

JMSPME_discrete_mass_plot(simulation_study)  ## find this function in auxiliary functions
