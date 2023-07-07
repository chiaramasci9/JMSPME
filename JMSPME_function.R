
JMSPME<-function(case,dataset,resp_var_name,ran_var_name=NULL,groups_name,iteration=30,iteration_for_support_reduction=20,iteration_fix=20,tollweight=0.05,toll=1,toll_ran_eff=10^-2,toll_fix_eff=10^-2) {
  library(nnet)
  caso=case
  IT=iteration #30 # total number of iterations
  IT1=iteration_for_support_reduction #20 #number of interation for the support reduction. From le iteration K1 to the end the support reduction is performed
  groups_var_name=groups_name
  #control parameters for the random effects distribution
  
  #distance between 2 points for make them collapse
  tollR=toll_ran_eff #tolerance for the iterative estimation of random effects support points
  
  #control parameters for the fixed effects estimate
  tollF=toll_fix_eff #tolerance for the iterative estimation of fixed effects
  itmax=iteration_fix #maximum number of iterations for the fixed effects estimate
  
  var_resp=dataset[,resp_var_name]
  groups=dataset[,groups_var_name]
  ind=c()
  ind[1]=which(names(dataset)==resp_var_name)
  ind[2]=which(names(dataset)==groups_var_name)
  if(!is.null(ran_var_name)){
    ind[3]=which(names(dataset)==ran_var_name)
  }
  fix_var_name=names(dataset)[-ind]
  n_fix=length(fix_var_name)
  
  if(caso=='int'){
    n_ran = 1
    random_param = 3
    if(is.null(ran_var_name)){
      ran_var=dataset[,fix_var_name[1]]
      fix_cov=as.data.frame(dataset[,fix_var_name[2:length(fix_var_name)]])
    }
    else{
      ran_var=dataset[,ran_var_name]
      fix_cov=as.data.frame(dataset[,fix_var_name])
    }
  }
  if(caso=='slope'){
    nrand = 1   #quanti sono i parametri random
    rand = 2
    ran_var=(dataset[,ran_var_name])
    fix_cov=as.data.frame(dataset[,fix_var_name])
  }
  if(caso=='inteslope'){
    nrand = 2   #quanti sono i parametri random
    rand = 12
    ran_var=dataset[,ran_var_name]
    fix_cov=as.data.frame(dataset[,fix_var_name])
  }
  
  groups=dataset[,groups_name]
  gruppi=levels(factor(groups))
  ind=c()
  for(i in 1:length(gruppi)){
    if(length(which(groups==gruppi[i]))>9) ind[length(ind)+1]=i #valutare il numero minimo di osservazioni come input
  }
  
  
  dati_ini=list()
  t_ini=list()
  n_rig_gruppoj=c()
  gruppi=gruppi[ind]
  for(j in 1:length(gruppi)){
    n_rig_gruppoj[j]=length(which(groups==gruppi[j]))
    k=which(groups==gruppi[j])
    for(i in 1:n_rig_gruppoj[j]){
      if(i==1) {
        dati_ini[[j]]=var_resp[k[i]]
        t_ini[[j]]=ran_var[k[i]]
      }
      else{
        dati_ini[[j]][i]=var_resp[k[i]]
        t_ini[[j]][i]=ran_var[k[i]]
      }
    }
  }
  
  n=length(gruppi)
  dim_max=max(n_rig_gruppoj)
  dati_n=matrix(NA,nrow=n,ncol=dim_max)
  t_n_2=matrix(NA,nrow=n,ncol=dim_max)
  for(i in 1:n){
    dati_n[i,1:n_rig_gruppoj[i]]=dati_ini[[i]]
    t_n_2[i,1:n_rig_gruppoj[i]]=t_ini[[i]]
  }
  n_fix=dim(as.data.frame(fix_cov))[2] #salvo il numero di covariate
  #creazione di una lista di liste
  t_fix=list()
  for(q in 1:n_fix){
    tn=list()
    for(j in 1:length(gruppi)){
      k=which(groups==gruppi[j])
      for(i in 1:n_rig_gruppoj[j]){
        if(i==1) {
          tn[[j]]=fix_cov[k[i],q]
        }
        else{
          tn[[j]][i]=fix_cov[k[i],q]
        }
      }
    }
    #a questo punto t_fix dovrebbe contenere in ogni elemento una lista che ha la j fix cov del gruppo i
    #ora trasformiamo in matrice
    t_n=matrix(NA,nrow=n,ncol=dim_max)
    for(i in 1:n){
      t_n[i,1:n_rig_gruppoj[i]]=tn[[i]]
    }
    
    t_fix[[q]]=t_n
    
  }
  dati=dati_n
  t=t_n_2
  t2=t_fix
  factor(dati[4,])
  
  dim(dati)
  dim(t)
  dim(t2)
  
  N=dim(dati)[1] 	#number of longitudinal curves
  nt=dim(dati)[2]
  
  
  
  K = length(levels(as.factor(as.vector(dati))))
  nparam = dim(dataf)[2]-1
  if(caso=='int'){
    nrand = 1   #quanti sono i parametri random
    rand = 1   #quali sono i parametri random (1 per intercetta, 2 per slope)
  }
  if(caso=='slope'){
    nrand = 1   #quanti sono i parametri random
    rand = 2   #quali sono i parametri random (1 per intercetta, 2 per slope)
  }
  if(caso=='inteslope'){
    nrand = 2   #quanti sono i parametri random
    rand = 12  #quali sono i parametri random (1 per intercetta, 2 per slope)
  }
  
  ####fixed and random effects range initialization
  
  
  # 1 is the reference category
  multi_mat <- matrix(NA,nrow=N,ncol=nparam*(K-1))
  converti_lista_in_matrice_tridimensionale <- function(lista_matrici) {
    num_matrici <- length(lista_matrici)
    righe <- nrow(lista_matrici[[1]])
    colonne <- ncol(lista_matrici[[1]])
    
    matrice_tridimensionale <- array(NA, dim = c(righe, colonne, num_matrici))
    
    for (i in 1:num_matrici) {
      matrice_tridimensionale[, , i] <- lista_matrici[[i]]
    }
    
    return(matrice_tridimensionale)
  }
  mat<-converti_lista_in_matrice_tridimensionale(t2)
  for(i in 1:N){
    yf = as.factor(dati[i,1:n_rig_gruppoj[i]])
    yf = relevel(yf, ref=1)
    mult = multinom(yf ~ t[i,1:n_rig_gruppoj[i]]+ mat[i,1:n_rig_gruppoj[i],], maxit = 1000)
    if(mult$convergence==0){
      multi_mat[i,1:nparam]=coefficients(mult)[1,]
      multi_mat[i,(nparam+1):(nparam*(K-1))]=coefficients(mult)[2,]
    }
    else{
      print("Warning: no convergence of the glm for the range initialization")
    }
  }
  # print(multi_mat)
  boxplot(multi_mat)
  range = list()
  
  ## alternativa 1: tengo tutti compresi i valori strani
  # for(k in 1:(K-1)){
  #   range[[k]]=list()
  #   for(p in 1:nparam){
  #     range[[k]][[p]] = c(range(multi_mat[,((k-1)*nparam + p)])[1]-1,
  #                         range(multi_mat[,((k-1)*nparam + p)])[2]+1)
  #   }
  # }
  
  
  for(k in 1:(K-1)){
    range[[k]]=list()
    for(p in 1:nparam){
      # range[[k]][[p]] = c(min(multi_mat[,((k-1)*nparam + p)], na.rm=T),
      #                     max(multi_mat[,((k-1)*nparam + p)], na.rm=T))
      
      # questo mi piaceva ma a volte ? troppo ampio
      ## alternativa 2: tengo solo fino ai valori dei baffi
      range[[k]][[p]] = c(quantile(multi_mat[,((k-1)*nparam + p)], 0.25, na.rm=T)-1.5*IQR(multi_mat[,((k-1)*nparam + p)], na.rm=T),
                          quantile(multi_mat[,((k-1)*nparam + p)], 0.75, na.rm=T)+1.5*IQR(multi_mat[,((k-1)*nparam + p)], na.rm=T))
    }
  }
  
  range
  
  ############ prova range a mano
  
  #range[[1]][[1]]=c(-4,1)
  #range[[1]][[2]]=c(-2,5)
  #range[[1]][[3]]=c(-5,-1)
  
  #range[[2]][[1]]=c(-4,2)
  #range[[2]][[2]]=c(-5,3)
  #range[[2]][[3]]=c(0,4)
  
  nknots=list()
  #### random effects initialization
  for(k in 1:(K-1)){
    if(N>20){
      nknots[[k]]=20
    }else{
      nknots[[k]]=N
    }
  }
  
  ## range initialization for nknots = N
  if(nknots[[1]]==N){
    knots = list()
    for(k in 1:(K-1)){
      if(rand==1) knots[[k]]= na.omit(multi_mat[,nparam*(k-1)+1])
      if(rand==2) knots[[k]]= na.omit(multi_mat[,nparam*(k-1)+2])
      if(rand==12) {
        knots[[k]]=list()
        knots[[k]][[1]]= na.omit(multi_mat[,nparam*(k-1)+1])
        knots[[k]][[2]]= na.omit(multi_mat[,nparam*(k-1)+2])
      }
      if(nrand==1) nknots[[k]]=length(knots[[k]])
      else if(nrand==2) nknots[[k]]=length(knots[[k]][[1]])
    }
  }
  
  
  ## range initialization for nknots < N
  if(nknots[[1]]<N){
    if(nrand==1){
      if(rand==1){
        knots=list()
        for(k in 1:(K-1)){
          knots[[k]]= (range[[k]][[1]][2] - range[[k]][[1]][1])*runif(nknots[[k]]) + range[[k]][[1]][1]
        }
      }
      else if(rand==2){
        knots=list()
        for(k in 1:(K-1)){
          knots[[k]]= (range[[k]][[2]][2] - range[[k]][[2]][1])*runif(nknots[[k]]) + range[[k]][[2]][1]
        }
      }
    }
    
    if(nrand==2){
      knots = list()
      
      for(k in 1:(K-1)){
        knots[[k]]=list()
        knots[[k]][[1]]= (range[[k]][[1]][2] - range[[k]][[1]][1])*runif(nknots[[k]]) + range[[k]][[1]][1]
        knots[[k]][[2]]= (range[[k]][[2]][2] - range[[k]][[2]][1])*runif(nknots[[k]]) + range[[k]][[2]][1]
      }
    }
  }
  
  
  #weights marginali
  weights_marg = list()
  for(k in 1:(K-1)){
    weights_marg[[k]]= rep(1/nknots[[k]],nknots[[k]])
  }
  
  weights_marg
  
  
  # weights congiunti
  dimw=nknots[[1]]
  for(k in 2:(K-1)){
    dimw= c(dimw,nknots[[k]])
  }
  dimw
  
  # partiamo con una distibuzione uniforme sulle masse
  someData <- rep(1/prod(dimw), prod(dimw))
  weights=array(someData,dimw)
  weights
  
  ## forse piuttosto farei la media (o la mediana?) dei coeff ottenuti in multi_mat, infatti:
  fix_param = list()
  if(nrand==1){
    for(k in 1:(K-1)){
      if(rand==1){
        fix_param[[k]] = colMeans(multi_mat[,(2+(nparam*(k-1))):(3+(nparam*(k-1)))],na.rm=T)
      }
      else if(rand==2){
        fix_param[[k]] = colMeans(multi_mat[,c(1+(nparam*(k-1)):3+(nparam*(k-1)))], na.rm=T)
      }
    }
  }
  if(nrand==2){
    for(k in 1:(K-1)){
      fix_param[[k]] = colMeans(as.matrix(multi_mat[,(3+(nparam*(k-1))) :(3+n_fix-1+(nparam*(k-1)))],ncol=n_fix), na.rm=T)
    }
  }
  
  
  
  #function definitions for the update
  #y is the vector of observations for the curve i
  
  #distance
  distE<-function(p1,p2){
    return(sqrt(sum((p1-p2)^2)))
  }
  
  ### pi_{ik} (? la likelihood di y_i) ci devo mettere la formula a pag 6 moltiplicata sui j,
  ### ma ho messo quella a pag 5
  
  
  ###############################################################################################
  
  
  
  
  ## facciamo la likelihood
  ## faccio direttamente l'esponenziale di eta, cos? posso fare il caso in cui non so a che gruppo appartengo con
  ## formula delle probabilit? totali
  
  likelihood <- function(index, knots, fix_param){
    # le categorie possibili sono da 2 a K
    # definisco la matrice di eta
    
    expeta = array(1, dim = c(N, dim_max , K))
    for(j in 1:n_rig_gruppoj[index]){
      for(kl in 2:K){
        
        if(rand==1){
          par1 = knots[[kl -1]]
          par2 = fix_param[[kl-1]][1]
          par3= fix_param[[kl-1]][2]
        }
        
        if(rand==2){
          par1 = fix_param[[kl-1]][1]
          par2 = knots[[kl -1]]
          par3= fix_param[[kl-1]][2]
        }
        
        if(rand==12){
          par1 = knots[[kl -1]][1]
          par2 = knots[[kl -1]][2]
          par3= fix_param[[kl-1]]
        }
        temp=0
        for(i in 1:n_fix){
          temp=temp+par3[i]*t2[[i]][index,j]
        }
        
        expeta[index, j, kl] = exp(par1 + par2*t[index,j] + temp)
        
      }
    }
    
    lik = 1
    for(j in 1:n_rig_gruppoj[index]){
      for(k in 1:K){ #prodotto su tutti i k, dove quella di riferimento ha eta=0
        lik = lik*ifelse(y[index,j]==k, (expeta[index,j,k])/sum(expeta[index,j,]),1)
      }
    }
    
    return(lik)
  }
  
  
  
  
  ### W_ij initialization
  compute_W <- function(){
    
    M = prod(unlist(nknots))
    #dimw = dim(weights) meglio mettere il numero di nknots
    dimW = c(N,unlist(nknots))
    someData <- rep(0, M)
    W = array(someData,dimW)
    W_temp = array(someData,dimW)
    
    
    for(m in 0:(M-1)){
      indexm = rep(0,(K-1))
      indexm[1]= floor(m/(prod(unlist(nknots)[-1])) )
      
      if (K==3)  indexm[2] = (m - indexm[1]*prod(unlist(nknots)[-1]))  # numeri interi per costruzione
      
      else if(K>3){
        for(k in 2:(K-1)){
          
          denomin = prod(unlist(nknots)[(k+1):(K-1)])
          
          iprec = indexm[1:(k-1)]
          
          Mprec = prod(unlist(nknots)[2:(K-1)])
          
          if(k>2){
            for(kk in 3: k){
              Mprec = c(Mprec,prod(unlist(nknots)[kk:(K-1)]))
            }
          }
          
          if(k < (K-1)) {
            indexm[k]= floor((m - sum(iprec*Mprec))/denomin)
          }
          else indexm[k]= floor(m - sum(iprec*Mprec))
        }
      }
      indexm = indexm +1 # aggiungo + 1 perch? il conto parte da 0
      
      knots_temp = list()
      if(nrand==1){
        for(k in 1:(K-1)){
          knots_temp[[k]]=knots[[k]][indexm[k]]
        }
      }
      if(nrand==2){
        for(k in 1:(K-1)){
          knots_temp[[k]]=c(knots[[k]][[1]][indexm[k]],knots[[k]][[2]][indexm[k]])
        }
      }
      
      # questo lo lascio cos? perch? non so come indicizzare automaticamente (regge da 3 a 5 categorie, ma basta aggiungerle)
      for(i in 1:N){
        if (K==3) W_temp[i,indexm[1],indexm[2]] <- weights[indexm[1],indexm[2]]*likelihood(i, knots_temp, fix_param) #
        if (K==4) W_temp[i,indexm[1],indexm[2],indexm[3]] <- weights[indexm[1],indexm[2],indexm[3]]*likelihood(i, knots_temp, fix_param) #
        if (K==5) W_temp[i,indexm[1],indexm[2],indexm[3],indexm[4]] <- weights[indexm[1],indexm[2],indexm[3],indexm[4]]*likelihood(i, knots_temp, fix_param) #
      }
    }
    
    
    for(i in 1:N){
      if(K==3) W_sum=sum(W_temp[i,,])
      else if(K==4) W_sum=sum(W_temp[i,,,])
      else if(K==5) W_sum=sum(W_temp[i,,,,])
      
      if(W_sum==0) W_sum = 1
      
      if(K==3) W[i,,]=W_temp[i,,]/W_sum
      else if(K==4) W[i,,,]=W_temp[i,,,]/W_sum
      else if(K==5) W[i,,,,]=W_temp[i,,,,]/W_sum
    }
    
    return(W)
    
  }
  
  
  W <- compute_W()
  W
  sum(W)
  
  
  compute_weights <- function(){
    weights = array(0, dim(W)[-1])
    
    M = (prod(unlist(nknots)))
    for(m in 0:(M-1)){
      indexm = rep(0,(K-1))
      indexm[1]= floor(m/(prod(unlist(nknots)[-1])) )
      
      if (K==3)  indexm[2] = (m - indexm[1]*prod(unlist(nknots)[-1]))  # numeri interi per costruzione
      
      else if(K>3){
        for(k in 2:(K-1)){
          
          denomin = prod(unlist(nknots)[(k+1):(K-1)])
          
          iprec = indexm[1:(k-1)]
          
          Mprec = prod(unlist(nknots)[2:(K-1)])
          
          if(k>2){
            for(kk in 3: k){
              Mprec = c(Mprec,prod(unlist(nknots)[kk:(K-1)]))
            }
          }
          
          if(k < (K-1)) {
            indexm[k]= floor((m - sum(iprec*Mprec))/denomin)
          }
          else indexm[k]= floor(m - sum(iprec*Mprec))
        }
      }
      indexm = indexm +1 # aggiungo + 1 perch? il conto parte da 0
      print(indexm)
      if(K==3) weights[indexm[1],indexm[2]] = sum(W[,indexm[1],indexm[2]])/N
      else if(K==4) weights[indexm[1],indexm[2],indexm[3]] = sum(W[,indexm[1],indexm[2],indexm[3]])/N
      else if(K==5) weights[indexm[1],indexm[2],indexm[3],indexm[4]] = sum(W[,indexm[1],indexm[2],indexm[3],indexm[4]])/N
      
    }
    return(weights)
  }
  
  weights = compute_weights()
  weights
  sum(weights)
  
  
  ## calcolo adesso le matrici W marginali per ogni k
  
  
  compute_W_marg <- function(){
    W_marg=list()
    for(k in 1:(K-1)){
      W_marg[[k]] = matrix(0, nrow=N, ncol= nknots[[k]])
      for(i in 1:N){
        for(nodo in 1:nknots[[k]]){
          
          if(k==1 & K==3) W_marg[[k]][i,nodo] = sum(W[i,nodo,])
          if(k==1 & K==4) W_marg[[k]][i,nodo] = sum(W[i,nodo,,])
          if(k==1 & K==5) W_marg[[k]][i,nodo] = sum(W[i,nodo,,,])
          
          if(k==2 & K==3) W_marg[[k]][i,nodo] = sum(W[i,,nodo])
          if(k==2 & K==4) W_marg[[k]][i,nodo] = sum(W[i,,nodo,])
          if(k==2 & K==5) W_marg[[k]][i,nodo] = sum(W[i,,nodo,,])
          
          if(k==3 & K==4) W_marg[[k]][i,nodo] = sum(W[i,,,nodo])
          if(k==3 & K==5) W_marg[[k]][i,nodo] = sum(W[i,,,nodo,])
          
          if(k==4 & K==5) W_marg[[k]][i,nodo] = sum(W[i,,,,nodo])
        }
        
      }
    }
    return(W_marg)
  }
  
  W_marg=compute_W_marg()
  W_marg
  sum(W_marg[[1]])
  
  compute_weights_marg <- function(){
    weights_marg = list()
    
    for(k in 1:(K-1)){
      weights_marg[[k]]= colSums(W_marg[[k]])/N
    }
    return(weights_marg)
  }
  
  weights_marg = compute_weights_marg()
  weights_marg
  sum(weights_marg[[1]])
  
  # eliminiamo i nodi che hanno pesi nulli
  for(k in 1:(K-1)){
    if(sum(colSums(W_marg[[k]])==0)>0){
      
      print("Warning: there is a zero-weight knots nella categoria:")
      print(k)
      Wnulli = which(colSums(W_marg[[k]])==0)
      
      W_marg[[k]]=W_marg[[k]][,-Wnulli]
      weights_temp_marg=weights_marg[[k]][-Wnulli]
      weights_marg[[k]]=weights_temp_marg/(sum(weights_temp_marg))
      nknots[[k]]=nknots[[k]]-length(Wnulli)
      
      if(nrand==1) knots[[k]]=knots[[k]][-Wnulli]
      else if(nrand==2) {
        knots[[k]][[1]]=knots[[k]][[1]][-Wnulli]
        knots[[k]][[2]]=knots[[k]][[2]][-Wnulli]
      }
      if(k==1 & K==3){
        W = W[,-Wnulli,]
        weights = weights[-Wnulli,]
      }
      if(k==1 & K==4){
        W = W[,-Wnulli,,]
        weights = weights[-Wnulli,,]
      }
      if(k==1 & K==5){
        W = W[,-Wnulli,,,]
        weights = weights[-Wnulli,,,]
      }
      
      if(k==2 & K==3){
        W = W[,,-Wnulli]
        weights = weights[,-Wnulli]
      }
      if(k==2 & K==4){
        W = W[,,-Wnulli,]
        weights = weights[,-Wnulli,]
      }
      if(k==2 & K==5){
        W = W[,,-Wnulli,,]
        weights = weights[,-Wnulli,,]
      }
      
      if(k==3 & K==4){
        W = W[,,,-Wnulli]
        weights = weights[,,-Wnulli]
      }
      if(k==3 & K==5){
        W = W[,,,-Wnulli,]
        weights = weights[,,-Wnulli]
      }
      
      if(k==4 & K==5){
        W = W[,,,,-Wnulli]
        weights = weights[,,,-Wnulli]
      }
      
    }
  }
  
  W
  dim(W)
  ####compute expectation random effects
  
  
  exp_c <- function(cr, W, nodo, categoria, knots, fix_param, weights){
    # categoria dovrebbe essere da 2 a K
    nknots_k = unlist(nknots)[-categoria]
    Mk = prod(nknots_k)
    valore = matrix(0, Mk, N)
    
    for(m in 0:(Mk-1)){
      if(K==3){
        indexm = rep(0,1)
        indexm[1] = m
      }
      else if(K>3){
        indexm = rep(0,(K-2))
        indexm[1]= floor(m/(prod(unlist(nknots_k)[-1])) )
        
        for(k in 2:(K-2)){
          
          denomin = prod(unlist(nknots_k)[(k+1):(K-1)])
          
          iprec = indexm[1:(k-1)]
          
          Mprec = prod(unlist(nknots_k)[2:(K-1)])
          
          if(k>2){
            for(kk in 3: k){
              Mprec = c(Mprec,prod(unlist(nknots_k)[kk:(K-1)]))
            }
          }
          
          if(k < (K-1)) {
            indexm[k]= floor((m - sum(iprec*Mprec))/denomin)
          }
          else indexm[k]= floor(m - sum(iprec*Mprec))
        }
      }
      indexm = indexm +1
      
      knots_temp = list()
      if(nrand==1){
        for(k in 1:(K-1)){
          if(k==categoria) knots_temp[[k]]= cr
          else{
            if(K==3) knots_temp[[k]]=knots[[k]][indexm[1]]
            else knots_temp[[k]]=knots[[k]][indexm[k]]
          }
        }
      }
      if(nrand==2){
        for(k in 1:(K-1)){
          if(k==categoria) knots_temp[[k]]=c(cr[1],cr[2])
          else{
            if(K==3) knots_temp[[k]]=c(knots[[k]][[1]][indexm[1]],knots[[k]][[2]][indexm[1]])
            else knots_temp[[k]]= c(knots[[k]][[1]][indexm[k]],knots[[k]][[2]][indexm[k]])
          }
        }
      }
      
      for(i in 1:N){
        if(categoria==1 & K==3) valore[m+1,i] = W[i,nodo, indexm[1]]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        else if(categoria==1 & K==4) valore[m+1,i] = W[i,nodo, indexm[1], indexm[2]]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        else if(categoria==1 & K==5) valore[m+1,i] = W[i,nodo, indexm[1], indexm[2], indexm[3]]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        
        else if(categoria==2 & K==3) valore[m+1,i] = W[i, indexm[1],nodo]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        else if(categoria==2 & K==4) valore[m+1,i] = W[i, indexm[1], nodo, indexm[2]]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        else if(categoria==2 & K==5) valore[m+1,i] = W[i,indexm[1], nodo, indexm[2], indexm[3]]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        
        else if(categoria==3 & K==4) valore[m+1,i] = W[i, indexm[1], indexm[2], nodo]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        else if(categoria==3 & K==5) valore[m+1,i] = W[i,indexm[1], indexm[2], nodo, indexm[3]]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
        
        else if(categoria==4 & K==5) valore[m+1,i] = W[i,indexm[1], indexm[2], indexm[3], nodo]*ifelse(likelihood(i,knots_temp, fix_param)==0,-exp(300),log(likelihood(i,knots_temp, fix_param)))
      }
    }
    
    # se dovesse servire
    # ifelse(likelihood(index, categoria, c[1], c[2], fix_param[[categoria-1]], knots, fix_param,weights)==0,
    # -exp(300),log(likelihood(index, categoria, c[1], c[2], fix_param[[categoria-1]], knots, fix_param, weights)))
    
    d = sum(valore)
    
    return(d)
  }
  
  
  
  ###to optimize random effects
  
  optim_ran <- function(){
    d = list()
    ran_hess = list()
    
    if(nrand==1){
      for(k in 1:(K-1)){
        d[[k]] <- rep(0,nknots[[k]])
        
        for(nodo in 1:nknots[[k]]){
          if(rand==1) range_opt = range[[k]][[1]]
          else if(rand==2) range_opt = range[[k]][[2]]
          categoria=k
          print(paste("categoria ", categoria))
          print(paste("nodo ", nodo))
          
          ottimn <- optim(knots_old[[k]][nodo], exp_c,exp_c, W=W, nodo=nodo, method="Brent",
                          categoria = k, lower = range_opt[1], upper = range_opt[2],
                          knots=knots_old, fix_param = fix_param, weights=weights_old,
                          control = list(fnscale = -1), hessian = T)
          d[[k]][nodo] <- ottimn$par
          ran_hess[[k]] <- ottimn$hessian
          #d[[k]][nodo] <- optimize(exp_c,range_opt, W=W, nodo = nodo,categoria = k, knots=knots_old, fix_param = fix_param, weights = weights_old, maximum = TRUE)$maximum
          
          if(exp_c(d[[k]][nodo],W, nodo, categoria=k, knots, fix_param, weights)==0) print("Warning: il massimo degli effetti random ? stato trovato in zero")
          
        }
      }
    }
    
    else if(nrand==2){
      for(k in 1:(K-1)){
        d[[k]]= list()
        d[[k]][[1]]= rep(0, nknots[[k]])
        d[[k]][[2]]= rep(0, nknots[[k]])
        range_opt1 = range[[k]][[1]]
        range_opt2 = range[[k]][[2]]
        for(nodo in 1:nknots[[k]]){
          ottimn = optim(c(knots_old[[k]][[1]][nodo],knots_old[[k]][[2]][nodo]), exp_c, W=W, nodo=nodo, method="Nelder-Mead",
                         categoria = k,
                         knots=knots_old, fix_param = fix_param, weights=weights_old,
                         control = list(fnscale = -1), hessian = T)
          ottim = ottimn$par
          ran_hess[[k]] = ottimn$hessian
          #ottim = optim(c(knots_old[[k]][[1]][nodo],knots_old[[k]][[2]][nodo]), exp_c, W=W, nodo=nodo, method="L-BFGS-B Nelder-Mead",
          #                lower=c(range_opt1[1],range_opt2[1]), upper =c(range_opt1[2],range_opt2[2]), categoria = k+1,
          #                knots=knots_old, fix_param = fix_param, weights=weights_old, control = list(fnscale = -1))$par
          
          if(exp_c(ottim,W, nodo, categoria=k, knots, fix_param, weights)==0) print("Warning: il massimo degli effetti random ? stato trovato in zero")
          
          d[[k]][[1]][nodo] <- ottim[1]
          d[[k]][[2]][nodo] <- ottim[2]
        }
      }
    }
    
    output=list()
    output[[1]]= d
    output[[2]]= ran_hess
    return(output)
  }
  
  ## voglio provare a disegnare la funzione da massimizzare nei random effects
  
  # for(k in 1:(K-1)){
  #   for(node in 1:nknots[[k]]){
  #     if(rand==1){
  #       xx = seq(from = range[[k]][[1]][1], to = range[[k]][[1]][2], by = 0.1)
  #       yy = rep(0, length(xx))
  #       for(i in 1:length(xx)){
  #         yy[i]=exp_c(xx[i], W, nodo=node, categoria=k+1, knots, fix_param, weights)
  #       }
  #       plot(xx, yy, main=paste('Nodo ', node,' della categoria ',k+1))
  #     }
  #
  #     if(rand==2){
  #       xx = seq(from = range[[k]][[2]][1], to = range[[k]][[2]][2], by = 0.1)
  #       yy = rep(0, length(xx))
  #       for(i in 1:length(xx)){
  #         yy[i]=exp_c(xx[i], W, nodo=node, categoria=k+1, knots, fix_param, weights)
  #       }
  #       plot(xx, yy, main=paste('Nodo ', node,' della categoria ',k+1))
  #     }
  #     if(rand==12){
  #       xx = seq(from = range[[k]][[1]][1], to = range[[k]][[1]][2], length = 100)
  #       yy = seq(from = range[[k]][[2]][1], to = range[[k]][[2]][2], length = 100)
  #       zz=rep(0, length(xx))
  #       for(i in 1:length(xx)){
  #         zz[i]=exp_c(c(xx[i], yy[i]), W, nodo=node, categoria=k+1, knots, fix_param, weights)
  #       }
  #       open3d()
  #       plot3d(0,0,0,type='n',xlim=range(xx),ylim=range(yy),
  #              zlim=range(zz))
  #       for(i in 1:length(zz)){
  #         #points3d(xx[i],yy[i],zz[i], col='black')
  #         points3d(xx[i],yy[i],zz[i], col='black', lwd=2)
  #       }
  #     }
  #   }
  # }
  
  
  
  
  #
  ####compute expectation beta
  
  
  exp_beta <- function(beta,knots,kbeta){
    M = prod(unlist(nknots))
    d <- rep(0,M)
    for(m in 0:(M-1)){
      
      indexm = rep(0,(K-1))
      indexm[1]= floor(m/(prod(unlist(nknots)[-1])) )
      
      if (K==3)  indexm[2] = (m - indexm[1]*prod(unlist(nknots)[-1]))  # numeri interi per costruzione
      
      else if(K>3){
        for(k in 2:(K-1)){
          
          denomin = prod(unlist(nknots)[(k+1):(K-1)])
          
          iprec = indexm[1:(k-1)]
          
          Mprec = prod(unlist(nknots)[2:(K-1)])
          
          if(k>2){
            for(kk in 3: k){
              Mprec = c(Mprec,prod(unlist(nknots)[kk:(K-1)]))
            }
          }
          
          if(k < (K-1)) {
            indexm[k]= floor((m - sum(iprec*Mprec))/denomin)
          }
          else indexm[k]= floor(m - sum(iprec*Mprec))
        }
      }
      indexm = indexm +1 # aggiungo + 1 perch? il conto parte da 0
      
      knots_temp = list()
      if(nrand==1){
        for(k in 1:(K-1)){
          knots_temp[[k]]=knots[[k]][indexm[k]]
          fix_param_temp = fix_param
          fix_param_temp[[kbeta]] = c(beta[1], beta[2])# qua modificare per i piu fix
        }
      }
      if(nrand==2){
        for(k in 1:(K-1)){
          knots_temp[[k]]=c(knots[[k]][[1]][indexm[k]],knots[[k]][[2]][indexm[k]])
          fix_param_temp = fix_param
          fix_param_temp[[kbeta]] = beta
        }
      }
      
      
      loglik = rep(0, N)
      for(i in 1:N){
        loglik[i] = ifelse(likelihood(i,knots_temp, fix_param_temp)==0,-exp(300),log(likelihood(i,knots_temp, fix_param_temp)))
      }
      
      # se dovesse servire
      # loglik[i] = ifelse(likelihood(i,categoria = kbeta+1, knots[[kbeta]][nodo],beta[1], beta[2], knots, fix_param,weights)==0,
      # -exp(300),log(likelihood(i,categoria = kbeta+1, knots[[kbeta]][nodo],beta[1], beta[2], knots, fix_param,weights)))
      
      if(K==3) W_vec = W[,indexm[1], indexm[2]]
      else if(K==4) W_vec = W[,indexm[1], indexm[2], indexm[3]]
      else if(K==5) W_vec = W[,indexm[1], indexm[2], indexm[3], indexm[4]]
      
      d[m] <- sum(W_vec*loglik)
    }
    return(sum(d))
  }
  
  
  
  ###optimize fixed variables
  
  optim_fixed <- function(){
    b=list()
    hess = list()
    if(nrand==1){
      for(k in 1:(K-1)){
        if(rand==1) {
          range_opt1 = range[[k]][[2]]
          range_opt2 = range[[k]][[3]]
        }
        if(rand==2) {
          range_opt1 = range[[k]][[1]]
          range_opt2 = range[[k]][[3]]
        }
        
        #b[[k]] <- optim(c(fix_param[[k]][1],fix_param[[k]][2]),exp_beta,knots=knots ,kbeta=k,
        #                 control = list(fnscale = -1))$par
        ottim = optim(c(fix_param[[k]]),exp_beta,knots=knots ,kbeta=k, method = "Nelder-Mead",
                      control = list(fnscale = -1), hessian=T)
        
        b[[k]] <- ottim$par
        hess[[k]] <- ottim$hessian
        #b[[k]] <- optim(c(fix_param[[k]][1],fix_param[[k]][2]),exp_beta,knots=knots ,kbeta=k, method = "L-BFGS-B",
        #                lower = c(range_opt1[1],range_opt2[1]), upper =c(range_opt1[2],range_opt2[2]), control = list(fnscale = -1))$par
        if(exp_beta(b[[k]],knots=knots, kbeta=k)==0) print("Warning: il massimo degli effetti fissi ? stato trovato in zero")
      }
    }
    if(nrand==2){
      for(k in 1:(K-1)){
        range_opt = range[[k]][[3]]
        if(n_fix==1){
          ottim <- optim(fix_param[[k]], exp_beta,knots=knots , lower = range_opt[1], upper=range_opt[2],
                         kbeta=k, method = "Brent", control = list(fnscale = -1), hessian=T)
        }
        if(n_fix>1) {
          ottim <- optim(fix_param[[k]], exp_beta,knots=knots , lower = range_opt[1], upper=range_opt[2],
                         kbeta=k, method = "Nelder-Mead", control = list(fnscale = -1), hessian=T)
        }
        b[[k]] <- ottim$par
        hess[[k]] <- ottim$hessian
        if(exp_beta(b[[k]],knots=knots, kbeta=k)==0) print("Warning: il massimo degli effetti fissi ? stato trovato in zero")
      }
    }
    output = list()
    output[[1]]=b
    output[[2]]=hess
    return(output)
  }
  
  
  
  IT=30 #30 # total number of iterations
  IT1=20 #20 #number of interation for the support reduction. From le iteration K1 to the end the support reduction is performed
  
  #control parameters for the random effects distribution
  tollweight=0.005
  
  toll=1 #distance between 2 points for make them collapse
  tollR=10^-2 #tolerance for the iterative estimation of random effects support points
  
  #control parameters for the fixed effects estimate
  tollF=10^-2 #tolerance for the iterative estimation of fixed effects
  itmax=5 #maximum number of iterations for the fixed effects estimate
  
  it=1
  conv1=0    #variabile ausiliaria che ci dice se ho raggiunto la convergenza rispetto alle tolleranze e quindi se posso applicare una support reduction
  conv2=0 #variabile ausiliaria che ci dice se ho cambiato il supporto e i rispettivi pesi oppure no
  colmean <- NULL
  
  knots_new=list()
  weights_new = weights
  weights_new_marg = list()
  weights_old_marg = list()
  
  
  while((conv1==0 | conv2==0) & (it < IT)){
    #support reduction of the discrete distribution support points based on the distance
    if(nrand==1){
      D=list()
      for(k in 1:(K-1)){
        D[[k]]=matrix(0,nknots[[k]],nknots[[k]])
        for(i in 1:(nknots[[k]]-1)){
          for(j in (i+1):nknots[[k]]){
            temp=distE(knots[[k]][i],knots[[k]][j])
            D[[k]][i,j]=temp
            D[[k]][j,i]=temp
          }
        }
        knots_new[[k]]=knots[[k]]
        weights_new_marg[[k]]=weights_marg[[k]]
        
        
        while(sum(D[[k]]<toll)>nknots[[k]]){
          
          id = which(D[[k]]==min(D[[k]][D[[k]]!=0]))[1]
          
          ####identifica riga e colonna del minimo valore D_ij
          
          if(id%%(dim(D[[k]])[1]) == 0) riga = floor(id/(dim(D[[k]])[1]))
          if(id%%(dim(D[[k]])[1]) != 0) riga = floor(id/(dim(D[[k]])[1])) + 1
          
          colonna = id -(riga-1)*(dim(D[[k]])[1])
          
          #knots_new[[k]][riga]=colMeans(rbind(knots_new[[k]][riga],knots_new[[k]][colonna]))
          knots_new[[k]][riga]= (knots_new[[k]][riga]*weights_new_marg[[k]][riga] + knots_new[[k]][colonna]*weights_new_marg[[k]][colonna])/(weights_new_marg[[k]][riga]+weights_new_marg[[k]][colonna])
          knots_new[[k]][colonna]=NA
          weights_new_marg[[k]][riga]=weights_new_marg[[k]][riga]+weights_new_marg[[k]][colonna]
          weights_new_marg[[k]][colonna]=NA
          
          nknots[[k]]=sum(!is.na(knots_new[[k]]))
          knots_new[[k]]=knots_new[[k]][!is.na(knots_new[[k]])]
          weights_new_marg[[k]]=weights_new_marg[[k]][!is.na(weights_new_marg[[k]])]
          
          # pesi congiunti
          if(k==1 & K==3){
            weights_new[riga,] = weights_new[riga,] + weights_new[colonna,]
            weights_new =weights_new[-colonna,]
          }
          if(k==1 & K==4){
            weights_new[riga,,] = weights_new[riga,,] + weights_new[colonna,,]
            weights_new =weights_new[-colonna,,]
          }
          if(k==1 & K==5){
            weights_new[riga,,,] = weights_new[riga,,,] + weights_new[colonna,,,]
            weights_new =weights_new[-colonna,,,]
          }
          
          if(k==2 & K==3){
            weights_new[,riga] = weights_new[,riga] + weights_new[,colonna]
            weights_new =weights_new[,-colonna]
          }
          if(k==2 & K==4){
            weights_new[,riga,] = weights_new[,riga,] + weights_new[,colonna,]
            weights_new =weights_new[,-colonna,]
          }
          
          if(k==2 & K==5){
            weights_new[,riga,,] = weights_new[,riga,,] + weights_new[,colonna,,]
            weights_new =weights_new[,-colonna,,]
          }
          
          if(k==3 & K==4){
            weights_new[,,riga] = weights_new[,,riga] + weights_new[,,colonna]
            weights_new =weights_new[,,-colonna]
          }
          
          if(k==3 & K==5){
            weights_new[,,riga,] = weights_new[,,riga,] + weights_new[,,colonna,]
            weights_new =weights_new[,,-colonna,]
          }
          
          if(k==4 & K==5){
            weights_new[,,,riga] = weights_new[,,,riga] + weights_new[,,,colonna]
            weights_new =weights_new[,,,-colonna]
          }
          
          
          D[[k]]=matrix(0,nknots[[k]],nknots[[k]])
          if(nknots[[k]]==1) {
            print("there are no different groups in confront")
            print(k)
            break
          }
          for(i in 1:(nknots[[k]]-1)){
            for(j in (i+1):nknots[[k]]){
              temp=distE(knots_new[[k]][i],knots_new[[k]][j])
              D[[k]][i,j]=temp
              D[[k]][j,i]=temp
            }
          }
          
        }
        
        weights_marg[[k]]=weights_new_marg[[k]]
        weights = weights_new
        knots[[k]]=knots_new[[k]]
        
      }
    }
    
    if(nrand==2){
      D=list()
      for(k in 1:(K-1)){
        D[[k]]=matrix(0,nknots[[k]],nknots[[k]])
        for(i in 1:(nknots[[k]]-1)){
          for(j in (i+1):nknots[[k]]){
            temp=distE(c(knots[[k]][[1]][i],knots[[k]][[2]][i]),c(knots[[k]][[1]][j],knots[[k]][[2]][j]))
            D[[k]][i,j]=temp
            D[[k]][j,i]=temp
          }
        }
        knots_new[[k]]=knots[[k]]
        weights_new_marg[[k]]=weights_marg[[k]]
        
        
        while(sum(D[[k]]<toll)>nknots[[k]]){
          
          id = which(D[[k]]==min(D[[k]][D[[k]]!=0]))[1]
          
          ####identifica riga e colonna del minimo valore D_ij
          
          if(id%%(dim(D[[k]])[1]) == 0) riga = floor(id/(dim(D[[k]])[1]))
          if(id%%(dim(D[[k]])[1]) != 0) riga = floor(id/(dim(D[[k]])[1])) + 1
          
          colonna = id -(riga-1)*(dim(D[[k]])[1])
          
          #knots_new[[k]][[1]][riga]=colMeans(rbind(knots_new[[k]][[1]][riga],knots_new[[k]][[1]][colonna]))
          #knots_new[[k]][[2]][riga]=colMeans(rbind(knots_new[[k]][[2]][riga],knots_new[[k]][[2]][colonna]))
          knots_new[[k]][[1]][riga]=(knots_new[[k]][[1]][riga]*weights_new_marg[[k]][riga] + knots_new[[k]][[1]][colonna]*weights_new_marg[[k]][colonna])/(weights_new_marg[[k]][riga]+weights_new_marg[[k]][colonna])
          knots_new[[k]][[2]][riga]=(knots_new[[k]][[2]][riga]*weights_new_marg[[k]][riga] + knots_new[[k]][[2]][colonna]*weights_new_marg[[k]][colonna])/(weights_new_marg[[k]][riga]+weights_new_marg[[k]][colonna])
          
          
          knots_new[[k]][[1]][colonna]=NA
          knots_new[[k]][[2]][colonna]=NA
          
          weights_new_marg[[k]][riga]=weights_new_marg[[k]][riga]+weights_new_marg[[k]][colonna]
          weights_new_marg[[k]][colonna]=NA
          
          nknots[[k]]=sum(!is.na(knots_new[[k]][[1]]))
          
          knots_new[[k]][[1]]=knots_new[[k]][[1]][!is.na(knots_new[[k]][[1]])]
          knots_new[[k]][[2]]=knots_new[[k]][[2]][!is.na(knots_new[[k]][[2]])]
          
          weights_new_marg[[k]]=weights_new_marg[[k]][!is.na(weights_new[[k]])]
          
          
          #nknots[[k]]=sum(!is.na(knots_new[[k]]))
          #knots_new[[k]]=knots_new[[k]][!is.na(knots_new[[k]])]
          #weights_new_marg[[k]]=weights_new_marg[[k]][!is.na(weights_new[[k]])]
          
          # pesi congiunti
          if(k==1 & K==3){
            weights_new[riga,] = weights_new[riga,] + weights_new[colonna,]
            weights_new =weights_new[-colonna,]
          }
          if(k==1 & K==4){
            weights_new[riga,,] = weights_new[riga,,] + weights_new[colonna,,]
            weights_new =weights_new[-colonna,,]
          }
          if(k==1 & K==5){
            weights_new[riga,,,] = weights_new[riga,,,] + weights_new[colonna,,,]
            weights_new =weights_new[-colonna,,,]
          }
          
          if(k==2 & K==3){
            weights_new[,riga] = weights_new[,riga] + weights_new[,colonna]
            weights_new =weights_new[,-colonna]
          }
          if(k==2 & K==4){
            weights_new[,riga,] = weights_new[,riga,] + weights_new[,colonna,]
            weights_new =weights_new[,-colonna,]
          }
          
          if(k==2 & K==5){
            weights_new[,riga,,] = weights_new[,riga,,] + weights_new[,colonna,,]
            weights_new =weights_new[,-colonna,,]
          }
          
          if(k==3 & K==4){
            weights_new[,,riga] = weights_new[,,riga] + weights_new[,,colonna]
            weights_new =weights_new[,,-colonna]
          }
          
          if(k==3 & K==5){
            weights_new[,,riga,] = weights_new[,,riga,] + weights_new[,,colonna,]
            weights_new =weights_new[,,-colonna,]
          }
          
          if(k==4 & K==5){
            weights_new[,,,riga] = weights_new[,,,riga] + weights_new[,,,colonna]
            weights_new =weights_new[,,,-colonna]
          }
          
          
          D[[k]]=matrix(0,nknots[[k]],nknots[[k]])
          if(nknots[[k]]==1) {
            print("there are no different groups in confront")
            print(k)
            break
          }
          for(i in 1:(nknots[[k]]-1)){
            for(j in (i+1):nknots[[k]]){
              temp=distE(c(knots_new[[k]][[1]][i],knots_new[[k]][[2]][i]),c(knots_new[[k]][[1]][j],knots_new[[k]][[2]][j]))
              D[[k]][i,j]=temp
              D[[k]][j,i]=temp
            }
          }
          
        }
        
        weights_marg[[k]]=weights_new_marg[[k]]
        weights=weights_new
        knots[[k]]=knots_new[[k]]
        
      }
    }
    
    
    
    #support reduction based on the weights
    if(conv1==1 | it>=IT1){
      nknots_old = nknots
      
      check = 0
      
      # group=list()
      # for(k in 1:(K-1)){
      # group[[k]]=rep(0,N)
      # for(i in 1:N){
      #   if (sum(W_marg[[k]][i,])==0) group[[k]][i]=NA
      #   else group[[k]][i]=which.max(W_marg[[k]][i,])
      # }
      
      group = matrix(0, N, (K-1) )
      for(i in 1:N){
        if(sum(W[i,,])!=0){
          if(K==3) group[i,] = which(W[i,,]==max(W[i,,]),arr.ind=T)
          else if(K==4) group[i,] = which(W[i,,,]==max(W[i,,,]),arr.ind=T)
          else if(K==5) group[i,] = which(W[i,,]==max(W[i,,,,]),arr.ind=T)
        }
      }
      
      for(k in 1:(K-1)){
        esisteNA = 0
        for(i in 1:nknots[[k]]){
          if( (weights_marg[[k]][i]<tollweight) & sum(group[,k]==i,na.rm=T)==0){
            print("Elimino dei nodi perch? hanno pesi bassissimi e non c'? nessuno")
            if(nrand==1) knots[[k]][i]=NA
            else if(nrand==2){
              knots[[k]][[1]][i]=NA
              knots[[k]][[2]][i]=NA
            }
            weights_marg[[k]][i]=NA
            esisteNA = esisteNA +1
          }
        }
        if(esisteNA>0){
          
          if(nrand==1) {
            nknots[[k]]=sum(!is.na(knots[[k]]))
            knots[[k]]=knots[[k]][!is.na(knots[[k]])]
          }
          else if(nrand==2){
            nknots[[k]]=sum(!is.na(knots[[k]][[1]]))
            knots[[k]][[1]]=knots[[k]][[1]][!is.na(knots[[k]][[1]])]
            knots[[k]][[2]]=knots[[k]][[2]][!is.na(knots[[k]][[2]])]
          }
          righe_na = which(is.na(weights_marg[[k]])) #qui c'? un problema perch? in teoria le ho gi? tolte le righe NA --> ho modificato  luglio 2021, vediamo
          print(paste("righe NA",righe_na))
          weights_marg[[k]]=weights_marg[[k]][!is.na(weights_marg[[k]])]
          
          
          if(K==3){
            if(k==1) weights = weights[-righe_na,]
            else if(k==2) weights = weights[,-righe_na]
          }
          else if(K==4){
            if(k==1) weights = weights[-righe_na,,]
            else if(k==2) weights = weights[,-righe_na,]
            else if(k==3) weights = weights[,,-righe_na]
          }
          else if(K==5){
            if(k==1) weights = weights[-righe_na,,,]
            else if(k==2) weights = weights[,-righe_na,,]
            else if(k==3) weights = weights[,,-righe_na,]
            else if(k==4) weights = weights[,,,-righe_na]
          }
          
          # ridistribuisco la massa persa
          weights_marg[[k]] = weights_marg[[k]]/(sum(weights_marg[[k]]))
          
          weights = weights/sum(weights)  # questo ? giusto pi? o meno, dovrei aggiornare weights marg degli altri k
          
          print(weights)
        }
        else check = check +1
      }
      
      if(check==(K-1)){
        conv2=1
        print("Ho raggiunto la convergenza 2")
      }
      print(paste('check',check))
      print('sono entrata in conv1 e i valori che ho alla fine di conv1 sono')
      print(paste('nknots', nknots))
      print(paste('knots', knots))
      print(paste('weights', weights))
      print(paste('weights_marg', weights_marg))
      print(paste('fixed', fix_param ))
    }
    # si dovrebbe ridistribuire la massa persa
    
    #computation of W_ij
    #W_old=W
    W=compute_W()
    dim(W)
    weights = compute_weights()
    
    W_marg=compute_W_marg()
    weights_marg = compute_weights_marg()
    
    
    # eliminiamo i nodi che hanno pesi nulli
    for(k in 1:(K-1)){
      if(sum(colSums(W_marg[[k]])==0)>0){
        
        print("Warning: there is a zero-weight knot nella categoria:")
        print(k)
        Wnulli = which(colSums(W_marg[[k]])==0)
        
        W_marg[[k]]=W_marg[[k]][,-Wnulli]
        weights_temp_marg=weights_marg[[k]][-Wnulli]
        weights_marg[[k]]=weights_temp_marg/(sum(weights_temp_marg))
        nknots[[k]]=nknots[[k]]-length(Wnulli)
        
        if(nrand==1) knots[[k]]=knots[[k]][-Wnulli]
        else if(nrand==2) {
          knots[[k]][[1]]=knots[[k]][[1]][-Wnulli]
          knots[[k]][[2]]=knots[[k]][[2]][-Wnulli]
        }
        if(k==1 & K==3){
          W = W[,-Wnulli,]
          weights = weights[-Wnulli,]
        }
        if(k==1 & K==4){
          W = W[,-Wnulli,,]
          weights = weights[-Wnulli,,]
        }
        if(k==1 & K==5){
          W = W[,-Wnulli,,,]
          weights = weights[-Wnulli,,,]
        }
        
        if(k==2 & K==3){
          W = W[,,-Wnulli]
          weights = weights[,-Wnulli]
        }
        if(k==2 & K==4){
          W = W[,,-Wnulli,]
          weights = weights[,-Wnulli,]
        }
        if(k==2 & K==5){
          W = W[,,-Wnulli,,]
          weights = weights[,-Wnulli,,]
        }
        
        if(k==3 & K==4){
          W = W[,,,-Wnulli]
          weights = weights[,,-Wnulli]
        }
        if(k==3 & K==5){
          W = W[,,,-Wnulli,]
          weights = weights[,,-Wnulli]
        }
        
        if(k==4 & K==5){
          W = W[,,,,-Wnulli]
          weights = weights[,,,-Wnulli]
        }
        
      }
    }
    
    # L'HO APPENA TOLTO PERCHE' SECONDO ME INUTILE
    # for(k in 1:(K-1)){
    # weights_old_marg[[k]]=weights_marg[[k]]
    # weights_old = weights
    # weights = compute_weights()
    # weights_marg = compute_weights_marg()
    # }
    
    weights_old = weights
    weights_old_marg=weights_marg
    #estimation of the random effects support points and fixed effects
    #iterative procedure
    itt=0
    
    fix_param_oldss = fix_param_old = fix_param
    knots_oldss = knots_old = knots
    
    
    #estimation of random effects support points
    knotsp <- optim_ran()
    knots <- knotsp[[1]]
    ran_hess <- knotsp[[2]]
    
    #fixed effects estimation
    fixed_temp=optim_fixed()
    fix_param =fixed_temp[[1]]
    fix_hess = fixed_temp[[2]]
    
    fix_param_diff=rep(0,(K-1)*(nparam-nrand))
    knots_diff = rep(0,sum(unlist(nknots)))
    
    for(k in 1:(K-1)){
      for(np in 1:(nparam-nrand)){
        fix_param_diff[(k-1)*(nparam-nrand) + np] = fix_param[[k]][np] - fix_param_old[[k]][np]
      }
    }
    
    indice=1
    for(k in 1:(K-1)){
      for(nnodi in 1:nknots[[k]]){
        if(nrand==1) knots_diff[indice] = knots[[k]][nnodi] - knots_old[[k]][nnodi]
        else if(nrand==2) knots_diff[indice] = (knots[[k]][[1]][nnodi] - knots_old[[k]][[1]][nnodi] + knots[[k]][[2]][nnodi] - knots_old[[k]][[2]][nnodi])/2
        
        indice =indice+1
      }
    }
    
    while( (max(abs(fix_param_diff))>tollF | max(abs(knots_diff))>tollR) & (itt<itmax) ){
      itt=itt+1
      fix_param_old = fix_param
      knots_old=knots
      
      #estimation of random effects support points
      
      knotsp <- optim_ran()
      knots <- knotsp[[1]]
      ran_hess <- knotsp[[2]]
      print(knots)
      #fixed effects estimation
      
      fixed_temp=optim_fixed()
      fix_param =fixed_temp[[1]]
      fix_hess = fixed_temp[[2]]
      print(fix_param)
      fix_param_diff=rep(0,(K-1)*(nparam-nrand))
      knots_diff = rep(0,sum(unlist(nknots)))
      
      for(k in 1:(K-1)){
        for(np in 1:(nparam-nrand)){
          fix_param_diff[(k-1)*(nparam-nrand) + np] = fix_param[[k]][np] - fix_param_old[[k]][np]
        }
      }
      
      indice=1
      for(k in 1:(K-1)){
        for(nnodi in 1:nknots[[k]]){
          if(nrand==1) knots_diff[indice] = knots[[k]][nnodi] - knots_old[[k]][nnodi]
          else if(nrand==2) knots_diff[indice] = (knots[[k]][[1]][nnodi] - knots_old[[k]][[1]][nnodi] + knots[[k]][[2]][nnodi] - knots_old[[k]][[2]][nnodi])/2
          
          indice =indice+1
        }
      }
      
      print("differenza tra knots")
      print(knots_diff)
      print("differenza tra fixed")
      print(fix_param_diff)
      print(itt)
    }
    
    
    #estimation of group, congiunto
    
    group = matrix(0, N, (K-1) )
    for(i in 1:N){
      if(sum(W[i,,])!=0){
        if(K==3) group[i,] = which(W[i,,]==max(W[i,,]),arr.ind=T)
        else if(K==4) group[i] = which(W[i,,,]==max(W[i,,,]),arr.ind=T)
        else if(K==5) group[i] = which(W[i,,]==max(W[i,,,,]),arr.ind=T)
      }
    }
    
    print("numero dell'iterazione")
    print(it)
    print(fix_param)
    print(knots)
    
    fix_param_diffss=rep(0,(K-1)*(nparam-nrand))
    knots_diffss = rep(0,sum(unlist(nknots)))
    
    for(k in 1:(K-1)){
      for(np in 1:(nparam-nrand)){
        fix_param_diffss[(k-1)*(nparam-nrand) + np] = fix_param[[k]][np] - fix_param_oldss[[k]][np]
      }
    }
    
    indice=1
    for(k in 1:(K-1)){
      for(nnodi in 1:nknots[[k]]){
        if(nrand==1) knots_diffss[indice] = knots[[k]][nnodi] - knots_oldss[[k]][nnodi]
        else if(nrand==2) knots_diffss[indice] = (knots[[k]][[1]][nnodi] - knots_oldss[[k]][[1]][nnodi] + knots[[k]][[2]][nnodi] - knots_oldss[[k]][[2]][nnodi])/2
        
        indice =indice+1
      }
    }
    
    if( max(abs(fix_param_diffss)<tollF) & (max(abs(knots_diffss)) < tollR) & conv1==0 ){
      conv1=1
      print("Prima convergenza raggiunta, ora inizio la riduzione dimensionale togliendo gli scarti")
    }
    
    it=it+1
  }
  
  
  
  
  
  return(list( Knots=knots, groups= group, weights=weights_marg, param=fix_param,case=caso))
  
}
