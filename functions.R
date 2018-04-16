#get log scale in x axis - no need
#Get everything in concentrations instead of number of molecules
#Do some thawing *
#Add the other kinetics - later. Let us not complicate things right now. 
library(gridExtra)
library(igraph)
library(RColorBrewer)
library(gplots)
library(rootSolve)
library(deSolve)



dt=0.01

updateRates<-function(reactions, metabolites, upperLimitConc){
  if(ncol(reactions)==5){
    reactions=cbind(reactions, matrix(0, nrow = nrow(reactions), ncol = 3))
    colnames(reactions)=c('reactant1', 'reactant2', 'product1', 'product2', 'rateConst', 'concReactant1', 'concReactant2', 'rate')
  }

  concentrations=metabolites[,ncol(metabolites)]
  reactions$concReactant1[reactions$reactant1>0]=concentrations[reactions$reactant1[reactions$reactant1>0]]
  reactions$concReactant1[reactions$reactant1==0]=rep(1, length(reactions$concReactant1[reactions$reactant1==0]))
  reactions$concReactant2[reactions$reactant2>0]=concentrations[reactions$reactant2[reactions$reactant2>0]]
  reactions$concReactant2[reactions$reactant2==0]=rep(1, length(reactions$concReactant2[reactions$reactant2==0]))
  reactions$rate=reactions$rateConst*reactions$concReactant1*reactions$concReactant2
  if (any(reactions$rate <0)){
    print('WAIT')
    print(metabolites)
    print(reactions)
  }
  return(reactions)
}

initNetwork<-function(NbSpecies, minrange, maxrange, upperLimitConc, lowerLimitConc, p=NULL, NbReactions=NULL){
  #print(p)
  metabolites<-as.data.frame(matrix(0, NbSpecies, 3))
  colnames(metabolites)<-c('name', 'mass', 'conc')
  metabolites$name<-1:NbSpecies
  metabolites$mass<-sample(minrange:maxrange, NbSpecies, replace = T)
  metabolites$conc=round((runif(length(metabolites$name), min = lowerLimitConc, max = upperLimitConc))*10000)/10000
  #changed
  metabolites<-metabolites[order(-metabolites$mass),]
  colnames<-c('reactant1', 'reactant2', 'product1', 'product2', 'rateConst', 'concReactant1', 'concReactant2', 'rate')
  reactions<-as.data.frame(matrix(0, 0, length(colnames)))
  colnames(reactions)<-colnames
  for(i in 1:(nrow(metabolites)-1)){
    #print(i)
    for(j in (i+1):(nrow(metabolites)-1)){
      if((metabolites$mass[i]-metabolites$mass[j]) %in% metabolites$mass[(j):nrow(metabolites)]){
        #print('---------------')
        #print(metabolites$name[i])
        #print(metabolites$name[j])
        #print((metabolites$mass[i]-metabolites$mass[j]))
        reactant2=metabolites$name[metabolites$mass==(metabolites$mass[i]-metabolites$mass[j])]
        #print(reactant2)
        #print('---------------')
        d=as.data.frame(matrix(0, length(reactant2), ncol(reactions)))
        colnames(d)<-colnames
        #reactions=rbind(reactions, matrix(0, length(reactant2), 4))
        #reactions$reactant1[(nrow(reactions)-length(reactant2)):nrow(reactions)]=metabolites$name[j]
        #reactions$reactant2[(nrow(reactions)-length(reactant2)):nrow(reactions)]=reactant2
        #reactions$product1[(nrow(reactions)-length(reactant2)):nrow(reactions)]=metabolites$name[i]
        d$reactant1[1:nrow(d)]=metabolites$name[j]
        d$reactant2[1:nrow(d)]=reactant2
        d$product1[1:nrow(d)]=metabolites$name[i]
        reactions=rbind(reactions, d)
        #do something
      }
    }
  }
  a=split(metabolites$name, metabolites$mass)
  pairs=t(as.data.frame(lapply(a[lapply(a, length)>1], function(x) combn(x,2))))
  #pairs
  if(nrow(pairs)>0){
    reactions[(nrow(reactions)+1):(nrow(reactions)+nrow(pairs)),]=cbind(pairs[,1], rep(0, nrow(pairs)), pairs[,2], rep(0, nrow(pairs)))
  } else{
    return(list(metabolites, reactions))
  }
  #print(nrow(reactions))
  reactions[(nrow(reactions)+1):(nrow(reactions)+nrow(reactions)),c(1,2,3,4,5,6,7,8)]=reactions[,c(3,4,1,2,5,6,7,8)]
  #View(reactions)
  reactions$rateConst=rep(1, nrow(reactions))
  #reactions$rateConst[(reactions$reactant1==0) | (reactions$reactant2==0)]=rep(1000, length(reactions$rateConst[(reactions$reactant1==0) | (reactions$reactant2==0)]))
  reactions$rateConst[(reactions$reactant1==0) | (reactions$reactant2==0)]=rep(1, length(reactions$rateConst[(reactions$reactant1==0) | (reactions$reactant2==0)]))
  if(!is.null(p)){
    reactions=reactions[unlist(lapply(runif(nrow(reactions), min = 0, max = 1), function(x) x>p)),]
  }
  else if(!is.null(NbReactions)){
    #print(nrow(reactions))
    #print(sample(1:nrow(reactions), NbReactions))
    reactions=reactions[sample(1:nrow(reactions), ifelse(NbReactions>nrow(reactions), nrow(reactions), NbReactions)),]
    #reactions=reactions[sample(1:nrow(reactions), NbReactions),]
    
  }
  
  #updating the reactions df with the concentrations of the reactants
  metabolites<-metabolites[order(metabolites$name),]
  reactions$concReactant1[reactions$reactant1>0]<-metabolites$conc[reactions$reactant1[reactions$reactant1>0]]
  reactions$concReactant2[reactions$reactant2>0]<-metabolites$conc[reactions$reactant2[reactions$reactant2>0]]
  reactions$concReactant1[reactions$reactant1==0]=rep(1, length(reactions$concReactant1[reactions$reactant1==0]))
  reactions$concReactant2[reactions$reactant2==0]=rep(1, length(reactions$concReactant2[reactions$reactant2==0]))
  return( list(metabolites, reactions))
}


"
findConstSpecies_old<-function(metabolites, reactions){
  #find species which are present only in reactants
  #find species which are present only in products
  #If the const species are still zero then we could add about 20% randomly chosen species.. but we can think about
  #that later
  SpeciesConst=NULL
  SpeciesConst=setdiff(metabolites$name, unique(c(reactions$product1, reactions$product2)))
  SpeciesConst=c(SpeciesConst, setdiff(metabolites$name, unique(c(reactions$reactant1, reactions$reactant2))))
  SpeciesConst=unique(SpeciesConst)
  return(SpeciesConst)
  
  #need to do some changes to this function
  #first make a network of metabolites
  #for each component, find the source and sink
  #if for any component, either the source or sink don't exist, then randomly assign the role to metabolites
  #return the array
}"


findConstSpecies<-function(metabolites, reactions, metaboliteNet=NULL){
  if(is.null(metaboliteNet)){
    metaboliteNet=getNetwork(metabolites = metabolites, reactions = reactions, print = F, display = F)$metaboliteNet
  }
  dg=decompose.graph(metaboliteNet)
  SpeciesConst=NULL
  for(i in 1:length(dg)){
    component=dg[[i]]
    if(length(vertex_attr(component)$name)==1){
      SpeciesConst=c(SpeciesConst, vertex_attr(component)$name)
      next
    }
    #find the indegree of each node of component
    #find the outdegree of each node of component
    #if none is zero for indegree then randoml choose some metbaoite
    #if none is zero for outdegree then reandomly choose some metabolite which is not the earlier one
    indegrees=degree(component, v=V(component), mode='in')
    sources=vertex_attr(component)$name[indegrees==0]
    outdegrees=degree(component, v=V(component), mode='out')
    sinks=vertex_attr(component)$name[outdegrees==0]
    if(length(sources)==0){
      #randomly choose one metabolite which is not in sinks
      sources=sample(setdiff(vertex_attr(component)$name, sinks), 1) 
    }
    if(length(sinks)==0){
      #randomly choose one metabolite which is not in sources
      sinks=sample(setdiff(vertex_attr(component)$name, sources), 1)
    }
    SpeciesConst=c(SpeciesConst, sources, sinks)
  }
  return(as.numeric(SpeciesConst))
}

plotConcProfile<-function(metabolites, log='', seed='NA', criticalPoints=NULL, sequence = NULL, maxConc=2){
  colnames<-colnames(metabolites)
  time=c(0, as.numeric(colnames[4:length(colnames)]))
  ymin=min(metabolites[,3:ncol(metabolites)])
  ymax=max(metabolites[,3:ncol(metabolites)])
  metabolites<-metabolites[order(metabolites$name),]
  
  #### Selecting colors
  n <- nrow(metabolites)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #par(mfrow=c(1,2))
  layout(matrix(c(1,1,1,1,2,1,1,1,1,2), 2, 5, byrow=TRUE))
  for(i in 1:nrow(metabolites)){
    if(i==1){
      plot(time, metabolites[i,3:ncol(metabolites)] ,col=col_vector[i], ylim = c(0, maxConc),  
           type = 'l', lwd=2, log=log,
           main = paste0('Concentration Profile \n Nb Species = ', as.character(NbSpecies), '\n Nb Reactions = ', as.character(nrow(reactions)), ' (p = ', as.character(p), ', seed = ', seed, ' )')
           ,ylab = 'Concentration', xlab = 'time'
           , cex.lab=1.5)
      lines(time, metabolites[i,3:ncol(metabolites)], col=col_vector[i])
      next
    }
    #points(time, metabolites[i,3:ncol(metabolites)], col=col_vector[i], pch=20)
    lines(time, metabolites[i,3:ncol(metabolites)], col=col_vector[i], lwd=2)
  }
  if(!(is.null(criticalPoints))){
    abline(v=criticalPoints, col='grey')
  }
  #if(!is.null(sequence)){
  #  sequence=c('Eqb',paste0('R', sequence))
  #  text(criticalPoints, rep(0, length(criticalPoints)), sequence)
  #}
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", ylab = '', xlab = '')
  legend( "top", legend=metabolites$name, col = col_vector[1:nrow(metabolites)], 
          lwd=2,title="Names", cex = 1.5)
}

old_select_perform_reaction<-function(reactions, metabolites, time, SpeciesConst, dt){
  #time=dt
  colnames<-colnames(metabolites)
  metabolites[,ncol(metabolites)+1]=metabolites[,ncol(metabolites)]
  colnames(metabolites)<-c(colnames, as.character(dt+time))
  
  for(index in 1:nrow(reactions)){
    reactant1=reactions$reactant1[index]
    reactant2=reactions$reactant2[index]
    product1=reactions$product1[index]
    product2=reactions$product2[index]
    
    dn=reactions$rate[index]*dt
    
    ######This is a hacky solution... Correct it
    if(reactant1>0){
      if((metabolites[reactant1, ncol(metabolites)]-dn)<=0){
        next
        #return(metabolites)
      }
      if(!(reactant1 %in% SpeciesConst)){
        metabolites[reactant1, ncol(metabolites)]=metabolites[reactant1, ncol(metabolites)]-dn
      }
    }
    if(reactant2>0){
      if((metabolites[reactant2, ncol(metabolites)]-dn)<=0){
        return(metabolites) 
      }
      if(!(reactant2 %in% SpeciesConst)){
        metabolites[reactant2, ncol(metabolites)]=metabolites[reactant2, ncol(metabolites)]-dn
      }
    }
    if(product1>0){
      if(!(product1 %in% SpeciesConst)){
        metabolites[product1, ncol(metabolites)]=metabolites[product1, ncol(metabolites)]+dn
      }
    }
    if(product2>0){
      if(!(product2 %in% SpeciesConst)){
        metabolites[product2, ncol(metabolites)]=metabolites[product2, ncol(metabolites)]+dn
      }
    }
  }
  
  return(metabolites)
}

################
# Matrix method of select_perform_reactions 
################

select_perform_reaction<-function(reactions, metabolites, time, SpeciesConst, dt=0.01, maxConc=2){
  colnames<-colnames(metabolites)
  metabolites[,(length(colnames)+1)]=metabolites[,length(colnames)]
  time=time+dt
  time=as.integer(time*10000)/10000
  
  colnames(metabolites)<-c(colnames, as.character(time))
  
  conc=metabolites[,ncol(metabolites)]
  S=matrixFromList(reactions, metabolites)
  #V=makeFluxMatrix(reactions)
  V=reactions$rate*dt
  dn=S%*%V
  dn[SpeciesConst]=0
  conc=conc+dn
  conc[conc<0]=0
  conc[conc<(10^(-7))]=0
  conc[conc>maxConc]=maxConc
  metabolites[,ncol(metabolites)]=conc
  
  return(metabolites) 
}

matrixFromList<-function(reactions, metabolites){
  rnmatrix=matrix(0, nrow(metabolites), nrow(reactions))
  for(i in 1:nrow(reactions)){
    reactants=reactions$reactant1[i]
    reactants=c(reactants, reactions$reactant2[i])
    products=reactions$product1[i]
    products=c(products, reactions$product2[i])
    for(j in 1:2){
      if(reactants[j]>0){
        rnmatrix[reactants[j], i]=rnmatrix[reactants[j], i]-1
      }
    }
    
    for(j in 1:2){
      if(products[j]>0){
        rnmatrix[products[j], i]=rnmatrix[products[j], i]+1
      }
    }
  }
  return(rnmatrix)
}


################
################

performReactions<-function(reactions, metabolites, maxTime=10000, intervals=100, maxCount=5000, dt=0.01,stopCriteria='Counts', SpeciesConst=NULL, upperLimitConc=5){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(metabolites, reactions)
  } 
  minCount=1000
  count=0
  if(ncol(metabolites)==3){
    time=0
  } else{
    time=as.numeric(colnames(metabolites)[ncol(metabolites)])
  }
  #while(cumTime<maxTime){
  #while(count<maxCount){
  stop=FALSE
  while(stop==FALSE){  
    time=time+dt
    #print(time)
    if(time%%2==0){
      print(time)
    }
    count=count+1
    reactions=updateRates(reactions, metabolites, upperLimitConc)
    metabolites=select_perform_reaction(reactions, metabolites, time, SpeciesConst, dt)
    if(count<minCount){
      #we dont want to check stopCriteria before enough number of timesteps have taken place
      next
    }
    if(stopCriteria=='time'){
      if(time>maxTime){
        stop=TRUE
      }
    }else if (stopCriteria=='Counts'){
      if(count>maxCount){
        stop=TRUE
      }
    }else if (stopCriteria=='Eqb'){
      #print(time)
      if(eqbExists(metabolites)){
        stop=TRUE
        print('eqb exists')
      }
    }
  }
  return(metabolites)
}

old_performReactions<-function(reactions, metabolites, maxTime=10000, intervals=100, maxCount=5000, dt=0.01, time=0, stopCriteria='Counts'){
  SpeciesConst=findConstSpecies(metabolites, reactions)
  count=0
  #while(cumTime<maxTime){
  #while(count<maxCount){
  stop=FALSE
  while(stop==FALSE){  
    time=time+dt
    print(time)
    count=count+1
    reactions=updateRates(reactions, metabolites, upperLimitConc)
    metabolites=select_perform_reaction(reactions, metabolites, time, SpeciesConst, dt)
    if(stopCriteria=='time'){
      if(time>Time){
        stop=TRUE
      }
    }else if (stopCriteria=='Counts'){
      if(count>maxCount){
        stop=TRUE
      }
    }else if (stopCriteria=='Eqb'){
      if(eqbExists(metabolites)){
        stop=TRUE
      }
    }
  }
  return(metabolites)
}

#plotFluxProfile<-function()

plotFluxProfile<-function(metabolites){
  metaboliteFlux=metabolites[,(1:(ncol(metabolites)-1))]
  for(i in 4:(ncol(metabolites)-1)){
    metaboliteFlux[,i]=metabolites[,(i+1)]-metabolites[,i]
  }
  colnames<-colnames(metaboliteFlux)
  time=c(0, as.numeric(colnames[4:length(colnames)]))
  ymin=min(metaboliteFlux[,3:ncol(metaboliteFlux)])
  ymax=max(metaboliteFlux[,3:ncol(metaboliteFlux)])
  metaboliteFlux<-metaboliteFlux[order(metaboliteFlux$mass),]
  
  #### Selecting colors
  n <- nrow(metaboliteFlux)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #par(mfrow=c(1,2))
  layout(matrix(c(1,1,1,1,2,1,1,1,1,2), 2, 5, byrow=TRUE))
  for(i in 1:nrow(metaboliteFlux)){
    if(i==1){
      plot(time, metaboliteFlux[i,3:ncol(metaboliteFlux)], col=col_vector[i], ylim = c(ymin, ymax),  
           type = 'l', lwd=2,
           main = paste0('Concentration Profile \n Nb Species = ', as.character(NbSpecies), '\n Nb Reactions = ', as.character(nrow(reactions)), ' (p = ', as.character(p), ' )')
           ,ylab = 'Concentration', xlab = 'time'
           , cex.lab=1.5)
      lines(time, metaboliteFlux[i,3:ncol(metaboliteFlux)], col=col_vector[i])
      next
    }
    #points(time, metaboliteFlux[i,3:ncol(metaboliteFlux)], col=col_vector[i], pch=20)
    lines(time, metaboliteFlux[i,3:ncol(metaboliteFlux)], col=col_vector[i], lwd=2)
  }
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", ylab = '', xlab = '')
  legend( "top", legend=metaboliteFlux$mass, col = col_vector[1:nrow(metaboliteFlux)], 
          lwd=2,title="Masses", cex = 1.5)
}

thaw <- function(sequence, metabolites, reactions, SpeciesConst=NULL){
  #uses thawSingleReaction()
  #uses eqbExists()
  maxConc=5
  colnames(metabolites) <- c('name', 'mass', 'conc', '0')
  #metabolites2=metabolites
  #reactions1=reactions
  #reactions=reactions1[1,]
  for(i in 1:length(sequence)){
    print(i)
    #metabolites<-thawSingleReaction(metabolites, reactions[sequence[i],], SpeciesConst)
    metabolites<-thawReactions(metabolites, reactions[sequence[1:i],], SpeciesConst)
    
    #modify selectPerformReactions such that it can run a reaction to its completion and can run an initiated reaction system.
    print('returned here as well')
    print(length(sequence))
  }
  plotConcProfile(metabolites)
  return(metabolites)
}

thawSingleReaction<-function(metabolites, reaction, SpeciesConst=NULL, n=1000, maxSteps=10000, maxConc=5000){
  #Run the reaction for n number of times
  #check for the last 0.01*n steps, if the slope is small
  #if slope is small(eqb achieved), return the metabolites
  #if slope is not small, run reaction for n steps more
  #what happens if the system starts oscillating, maybe just print that something is wrong with 'this' reaction and continue with some other reaction. 
  print('welcome')
  reactant1=reaction$reactant1[1]
  reactant2=reaction$reactant2[1]
  product1=reaction$product1[1]
  product2=reaction$product2[1]
  
  time=as.numeric(colnames(metabolites)[ncol(metabolites)])
  i=1
  print(time)
  while(i<=n & i<maxSteps){
    #print(i)
    if(i==n){
      #check if eqb make a new function
      print('i reaches n')
      if(eqbExists(metabolites)){
        print('reaches eqb')
        break
      }
      else{
        n=n+500
        print('proceeds')
      }
    }
    #perform single reaction
    #go to next reaction if eqb reached.
    dn=reaction$rate[1]*dt
    time=time+dt
    time=as.integer(time*10000)/10000
    colnames=colnames(metabolites)
    metabolites=cbind(metabolites, metabolites[,ncol(metabolites)])
    colnames(metabolites)=c(colnames, as.character(time))
    flag=0
    if(reactant1>0 & !(reactant1 %in% SpeciesConst)){
      if((metabolites[reactant1, ncol(metabolites)]-dn)<0){
        #print('reaches min')
        flag=1
      }
    }
    if(reactant2>0 & !(reactant2 %in% SpeciesConst)){
      if((metabolites[reactant2, ncol(metabolites)]-dn)<0){
        #print('reaches min')
        flag=1
      }
    }
    if(product1>0 & !(product1 %in% SpeciesConst)){
      if((metabolites[product1, ncol(metabolites)]+dn)>maxConc){
        #print('reaches max')
        flag=1
      }
    }
    if(product2>0 & !(product2 %in% SpeciesConst)){
      if((metabolites[product2, ncol(metabolites)]+dn)>maxConc){
        #print('reaches max')
        flag=1
      }
    }
    
    if(flag==1){
      i=i+1
      next
    }
    
    if(reactant1>0){
      if(!(reactant1 %in% SpeciesConst)){
        metabolites[reactant1, ncol(metabolites)]=metabolites[reactant1, ncol(metabolites)]-dn
      }
    }
    if(reactant2>0){
      if(!(reactant2 %in% SpeciesConst)){
        metabolites[reactant2, ncol(metabolites)]=metabolites[reactant2, ncol(metabolites)]-dn
      }
    }
    if(product1>0){
      if(!(product1 %in% SpeciesConst)){
        metabolites[product1, ncol(metabolites)]=metabolites[product1, ncol(metabolites)]+dn
      }
    }
    if(product2>0){
      if(!(product2 %in% SpeciesConst)){
        metabolites[product2, ncol(metabolites)]=metabolites[product2, ncol(metabolites)]+dn
      }
    }
    
    i=i+1
    
    if(i%%1000==0){
      print(i)
      print(reaction)
    }
  }
  print(i)
  print('return')
  return(metabolites)
}

eqbExists<-function(metabolites, testSteps=1000, maxConc=NULL){
  metabolites=data.frame(metabolites, stringsAsFactors = F)
  if(ncol(metabolites)<testSteps){
    return(FALSE)
  }
  if(is.null(maxConc)){
    #maxConc=max(metabolites[,4:ncol(metabolites)])
    maxConc=5
    #print('hey')
  }
  #check if the last testSteps conc of each metabolite doesn't change by 10% of maxConc
  sdConc=unlist(lapply(metabolites$name, function(x) sd(metabolites[x,(ncol(metabolites)-testSteps):ncol(metabolites)])))
  #print('bye')
  if(any(sdConc>0.001*maxConc)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

thawReactions <- function(metabolites, reactions, SpeciesConst=NULL, dt=0.01,
                          n=1000, maxSteps=10000, maxConc=5000){
  
  ############ 
  #Right now this code runs for just one reaction or something like that
  
  print(reactions)
  #while i<=n and i<max{
  #run all reactions (j : 1:nrow(reactions))
  #if i==n check for the last 0.01*n steps, if the slope is small
  #if slope is small(eqb achieved), return the metabolites
  #if slope is not small, run reaction for n steps more
  #}
  i=1
  while(i<=n & i<max){
    
    if(i==n){
      #check if eqb make a new function
      print('i reaches n')
      if(eqbExists(metabolites)){
        print('reaches eqb')
        break
      }
      else{
        n=n+500
        print('proceeds')
      }
    }
    
    time=as.numeric(colnames(metabolites)[ncol(metabolites)])
    #print(time)
    time=time+dt
    time=as.integer(time*10000)/10000
    colnames=colnames(metabolites)
    metabolites=cbind(metabolites, metabolites[,ncol(metabolites)])
    colnames(metabolites)=c(colnames, as.character(time))
    
    for(index in 1:nrow(reactions)){
      reactant1=reactions$reactant1[index]
      reactant2=reactions$reactant2[index]
      product1=reactions$product1[index]
      product2=reactions$product2[index]
      
      dn=reactions$rate[index]*dt
      
      flag=0
      if(reactant1>0 & !(reactant1 %in% SpeciesConst)){
        if((metabolites[reactant1, ncol(metabolites)]-dn)<0){
          #print('reaches min')
          flag=1
        }
      }
      if(reactant2>0 & !(reactant2 %in% SpeciesConst)){
        if((metabolites[reactant2, ncol(metabolites)]-dn)<0){
          #print('reaches min')
          flag=1
        }
      }
      if(product1>0 & !(product1 %in% SpeciesConst)){
        if((metabolites[product1, ncol(metabolites)]+dn)>maxConc){
          #print('reaches max')
          flag=1
        }
      }
      if(product2>0 & !(product2 %in% SpeciesConst)){
        if((metabolites[product2, ncol(metabolites)]+dn)>maxConc){
          #print('reaches max')
          flag=1
        }
      }
      
      if(flag==1){
        next
      }
      
      if(reactant1>0){
        if(!(reactant1 %in% SpeciesConst)){
          metabolites[reactant1, ncol(metabolites)]=metabolites[reactant1, ncol(metabolites)]-dn
        }
      }
      if(reactant2>0){
        if(!(reactant2 %in% SpeciesConst)){
          metabolites[reactant2, ncol(metabolites)]=metabolites[reactant2, ncol(metabolites)]-dn
        }
      }
      if(product1>0){
        if(!(product1 %in% SpeciesConst)){
          metabolites[product1, ncol(metabolites)]=metabolites[product1, ncol(metabolites)]+dn
        }
      }
      if(product2>0){
        if(!(product2 %in% SpeciesConst)){
          metabolites[product2, ncol(metabolites)]=metabolites[product2, ncol(metabolites)]+dn
        }
      }
      
    }
    
    i=i+1
    
    if(i%%1000==0){
      print(i)
      print(reactions)
    }
    
  }
  
  return(metabolites)
}

"
thaw_oneByone<-function(metabolites, reactions, SpeciesConst=NULL){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies_old(reactions, metabolites)
  }
  active_rn=NULL
  criticalPoints=NULL
  criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  while(length(active_rn)<nrow(reactions)){
    #switching on rn randomly one by one
    active_rn=c(active_rn, sample(setdiff((1:nrow(reactions)), active_rn), 1))
    metabolites=performReactions(reactions[active_rn,], metabolites, SpeciesConst = SpeciesConst, stopCriteria = 'Eqb')
    criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  }
  return(list(metabolites, criticalPoints))
}
"
thaw_oneByone_runsteady<-function(metabolites, eqbProfile=NULL, reactions, SpeciesConst=NULL, maxConc=2){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(reactions, metabolites)
  }
  active_rn=NULL
  if(is.null(eqbProfile)){
    eqbProfile=metabolites[, ncol(metabolites)]
  }
  profiles<-matrix(eqbProfile, nrow=1)
  while(length(active_rn)<nrow(reactions)){
    active_rn=c(active_rn, sample(setdiff((1:nrow(reactions)), active_rn), 1))
    S=matrixFromList(reactions = reactions[active_rn,], metabolites = metabolites)
    yini<-profiles[nrow(profiles),]
    out<-runsteady(time=c(0,1e5), y=yini, func=fluxCalc, 
                   parms=list(S, SpeciesConst, maxConc, reactions$rateConst[active_rn]))
    profiles<-rbind(profiles, out$y)
  }
  return(profiles[nrow(profiles),])
}

thaw_inBulk<-function(metabolites, reactions, SpeciesConst=NULL, fraction=0.5){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(reactions, metabolites)
  }
  active_rn=NULL
  criticalPoints=NULL
  criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  for(i in 1:2){
    #switching on rn randomly one by one
    if(i==1){
      active_rn=c(active_rn, sample(setdiff((1:nrow(reactions)), active_rn), as.integer(fraction*nrow(reactions)), replace = FALSE))
    }
    else{
      active_rn=1:nrow(reactions)
    }
    metabolites=performReactions(reactions[active_rn,], metabolites, SpeciesConst = SpeciesConst, stopCriteria = 'Eqb')
    criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  }
  return(list(metabolites, criticalPoints))
}


"fluxCalc<-function(t, y, parms){
  #parms[[1]]<-S -> this is constant
  #parms[[2]]<-SpeciesConst
  #parms[[3]]<-maxConc
  #parms[[4]]<-rateConsts
  #conc<-parms[[1]]
  S<-parms[[1]]
  SpeciesConst<-parms[[2]]
  maxConc<-parms[[3]]
  rateConsts<-parms[[4]]
  rate=unlist(lapply(1:ncol(S), function(x) prod(y[S[,x]<0])))
  rate=rate*rateConsts
  dy=S%*%matrix(rate, nrow = length(rate), ncol = 1)
  inds1<-which((y+dy)<0)
  inds2<-which((y+dy)>maxConc)
  #inds1<-which(((y+dy)<1e-7)&(dy<0)) 
  #inds2<-which(((y+dy)>maxConc)&(dy>0))
  if(length(inds1)>0){
    dy[inds1]=-(y[inds1]) #this will convert negative y into 0
  }
  if(length(inds2)>0){
    dy[inds2]=(rep(maxConc, length(inds2))-y[inds2]) #this will bound the conc by an upperbound
  }
  dy[SpeciesConst]=0
  #simplify matrix to vector
  list(as.vector(dy))
}"

fluxCalc<-function(t, y, parms){
  #merge this into the fluxCalc function once we have compared.
  S<-parms[[1]]
  SpeciesConst<-parms[[2]]
  maxConc<-parms[[3]]
  rateConsts<-parms[[4]]
  minConc=1e-4
  #print(rateConsts)
  #if(length(parms)==4){
  #  Kms=matrix(rep(1, nrow(S)*ncol(S)), nrow(S), ncol(S)) # just setting all Kms to 1 for now
  #} else {
  #  Kms<-parms[[5]]
  #}
  #rate=unlist(lapply(1:ncol(S), function(x) prod(y[S[,x]<0])))
  #print(rateConsts)
  rate=convenience(S = S, y = y, kcats = rateConsts)
  #rate=rate*rateConsts
  dy=sign(S)%*%matrix(rate, nrow = length(rate), ncol = 1) #not to be multiplied by the coefficient , just by the sign of the element
  inds1<-which((y+dy)<minConc)
  inds2<-which((y+dy)>maxConc)
  #inds1<-which(((y+dy)<1e-7)&(dy<0)) 
  #inds2<-which(((y+dy)>maxConc)&(dy>0))
  if(length(inds1)>0){
    dy[inds1]=minConc-(y[inds1]) #this will convert negative y into 0
  }
  if(length(inds2)>0){
    dy[inds2]=(rep(maxConc, length(inds2))-y[inds2]) #this will bound the conc by an upperbound
  }
  dy[SpeciesConst]=0
  #simplify matrix to vector
  #print(as.vector(dy))
  list(as.vector(dy))
}

convenience<-function(S, y, kcats=NULL, Kms=NULL, E_conc=NULL){
  #E_conc is the concentration of the enzyme which we are assuming to be 1
  #kcats are assumed to be kcat+
  #right now this function just assumes irreversible 
  species=1:nrow(S)
  if(is.null(E_conc)){
    E_conc=1
  }
  if(is.null(kcats)){
    kcats=rep(1, ncol(S))
  }
  if(is.null(Kms)){
    Kms=S/S
  }
  rates=NULL
  for(i in 1:ncol(S)){
    reactants=species[S[,i]<0]
    products=species[S[,i]>0]
    yreactants_=y[reactants]/Kms[reactants,i]
    yproducts_=y[products]/Kms[products,i]
    coeffs=abs(S[reactants,i])
    numerator=E_conc*kcats[i]*prod(yreactants_^coeffs)
    denominator=prod(unlist(lapply(1:length(reactants), function(x) sum(yreactants_[x]^(0:coeffs[x]))))) 
    rates=c(rates, numerator/denominator)
  }
  return(rates)
}

NbBranches<-function(metabolites,reactions){
  #for each metabolite
  ##if the metabolite is in the product
  ##x is the number of reactions in which it is a reactant
  ##x-1 is the number of branching in this metabolite, add this to the total number of branches
  branches=0
  for(i in 1:nrow(metabolites)){
    if(i %in% union(reactions$product1, reactions$product2)){
      inds<-union(which(reactions$reactant1 %in% i), which(reactions$reactant2 %in% i))
      inds<-setdiff(inds, union(which(reactions$product1 %in% i), which(reactions$product2 %in% i)))
      if(length(inds)>1){
        branches = branches + length(inds) - 1
      }
    }
  }
  return(branches)
}

NbBranchingMetabolites<-function(metabolites,reactions){
  #for each metabolite
  ##if the metabolite is in the product
  ##x is the number of reactions in which it is a reactant
  ##if x>1 then add that metabolite to the list of branching metabolites. 
  branchingMet=0
  for(i in 1:nrow(metabolites)){
    if(i %in% union(reactions$product1, reactions$product2)){
      inds<-union(which(reactions$reactant1 %in% i), which(reactions$reactant2 %in% i))
      inds<-setdiff(inds, union(which(reactions$product1 %in% i), which(reactions$product2 %in% i)))
      if(length(inds)>1){
        branchingMet = branchingMet + 1
      }
    }
  }
  return(branchingMet)
}


getNetwork<-function(metabolites, reactions, print=F, resultDir=NULL, filename='reaction_network.pdf', display = F, SpeciesConst=NULL){
  rownames(reactions)=NULL
  bipartite<-NULL
  reactionNet<-NULL
  metaboliteNet<-NULL
  df=matrix(0, 0, 3)
  for(i in 1:nrow(reactions)){
    reactants=reactions[i, c(1,2)]
    products=reactions[i,c(3,4)]
    reactants=reactants[reactants>0]
    products=products[products>0]
    df<-rbind(df, cbind(expand.grid(reactants, products), i))
  }
  df=cbind(df, paste0('R',as.character(df[,3])))
  df=data.frame(df, stringsAsFactors = F)
  colnames(df)=c('V1', 'V2', 'R', 'R_name')
  
  #######METABOLITE NETWORK
  #######
  metaboliteNet<-make_graph(as.character(as.vector(t(df[,c(1,2)]))), 
                            isolates = as.character(setdiff(1:nrow(metabolites),unique(unlist(as.vector(reactions[, c(1,2,3,4)]))))) ,
                            directed = T)
  #plot(metaboliteNet)
  #######
  
  #######BIPARTITE NETWORK
  #########
  #bipartite=make_bipartite_graph(rep(c(0,1), length=length(edges)), as.character(edges))
  g <- graph.empty()
  g <- add.vertices(g,nv=nrow(reactions),attr=list(name=paste0('R',1:nrow(reactions)),
                                                   type=rep(1,nrow(reactions))))
  g <- add.vertices(g,nv=nrow(metabolites),attr=list(name=as.character(1:nrow(metabolites)),
                                                     type=rep(0,nrow(metabolites))))
  edges=as.character(c(as.vector(t(unique(df[,c(1,4)]))), as.vector(t(unique(df[,c(4,2)])))))
  edges=gsub(' ', '', edges) #hacky solution to remove weird white spaces being added in the vertex names
  # we need to turn edgeList into a vector (and using names instead of indexes)
  #edgeListVec <- as.vector(t(as.matrix(data.frame(S1=paste0('A',edgeList$S1),
  #                                                S2=paste0('B',edgeList$S2)))))
  g<- add.edges(g,edges)
  
  #plot.igraph(g, 
  #            vertex.color=c("orange","green")[as.numeric(V(g)$type)+1],
  #            vertex.shape = c('circle', 'square')[as.numeric(V(g)$type)+1])
  
  bipartite<-g
  ############
  
  #######REACTION NETWORK
  #########
  #if products of one reaction is a reactant in another reaction then make a directed edge from first network to second network
  #for each met, get the reactions in which it is a reactant
  #              get the reaction in which it is a product
  #              if the length of both >0 hten expand the two arrays
  pairs=NULL
  for(i in 1:nrow(metabolites)){
    reactions1=unique(df$R_name[df$V1==i])
    reactions2=unique(df$R_name[df$V2==i])
    pairs=rbind(pairs, expand.grid(reactions2, reactions1))
  }
  pairs=unique(pairs)
  reactionNet=make_graph(as.character(as.vector(t(pairs))), 
                         isolates = setdiff(df$R_name, unique(as.vector(t(pairs)))) ,
                         directed = T)
  ########
  
  ########PLOTTING
  #######
  if(print==T){
    SpeciesConst=SpeciesConst
    #print(SpeciesConst)
    print(paste0(resultDir, filename))
    pdf(paste0(resultDir, filename))
    layout=matrix(c(0,0,1,1,1,1,0,0), 4, 8, byrow = T)
    layout=rbind(layout, matrix(c(2,2,2,2,3,3,3,3), 4, 8, byrow=T))
    layout(layout)
    g<-bipartite
    colors<-c("orange","green")[as.numeric(V(g)$type)+1]
    colors[V(g)$name %in% as.character(SpeciesConst)]='red'
    plot.igraph(g, 
                vertex.color=colors,
                vertex.shape = c('circle', 'square')[as.numeric(V(g)$type)+1],
                #layout=layout,
                edge.arrow.size=.4, edge.curved=.3,
                main = paste0('Nb metabolites = ', nrow(metabolites), '\nNb reactions = ', nrow(reactions),
                              '\nseed = ', seed), cex.main=3)
    layout=layout_with_fr(metaboliteNet)*10
    colors=rep('orange', length(V(metaboliteNet)))
    colors[V(metaboliteNet)$name %in% as.character(SpeciesConst)]='red'
    plot(metaboliteNet, #layout=layout,
         vertex.color=colors,
         edge.arrow.size=.4, edge.curved=.3)
    layout=layout_with_fr(reactionNet)*10
    plot(reactionNet, vertex.shape='square', vertex.color='green', #layout=layout,
         edge.arrow.size=.6, edge.curved=.3)
    plot.new()
    grid.table(reactions[1:4]) #nextpage
    
    dev.off()
  }
  if(display==T){    
    SpeciesConst=SpeciesConst
    layout=matrix(c(0,0,1,1,1,1,0,0), 4, 8, byrow = T)
    layout=rbind(layout, matrix(c(2,2,2,2,3,3,3,3), 4, 8, byrow=T))
    layout(layout)
    
    #par(mfrow=c(2,2))
    g<-bipartite
    colors<-c("orange","green")[as.numeric(V(g)$type)+1]
    colors[V(g)$name %in% as.character(SpeciesConst)]='red'
    plot.igraph(g, 
                vertex.color=colors,
                vertex.shape = c('circle', 'square')[as.numeric(V(g)$type)+1],
                #layout=layout,
                edge.arrow.size=.4, edge.curved=.3,
                main = paste0('Nb metabolites = ', nrow(metabolites), '\nNb reactions = ', nrow(reactions),
                              '\nseed = ', seed), cex.main=3)
    layout=layout_with_fr(metaboliteNet)*10
    colors=rep('orange', length(V(metaboliteNet)))
    colors[V(metaboliteNet)$name %in% as.character(SpeciesConst)]='red'
    plot(metaboliteNet, #layout=layout, 
         vertex.color=colors,
         edge.arrow.size=.4, edge.curved=.3)
    layout=layout_with_fr(reactionNet)*10
    plot(reactionNet, vertex.shape='square', vertex.color='green', #layout=layout,
         edge.arrow.size=.6, edge.curved=.3)
  }
  #######
  
  return(list(bipartite =  bipartite, reactionNet = reactionNet, metaboliteNet = metaboliteNet))
}

'
thaw_oneByone_ode<-function(metabolites, reactions, SpeciesConst=NULL, maxConc=2, print = F){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(reactions, metabolites)
  }
  active_rn=NULL
  criticalPoints=NULL
  criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  reactions_full=reactions
  df=data.frame(matrix(0, nrow(metabolites), (nrow(reactions)+1)))
  colnames(df)=c("Eqb", c("R", 1:nrow(Reactions)))
  while(length(active_rn)<nrow(reactions_full)){
    #switching on rn randomly one by one
    #print(setdiff((1:nrow(reactions_full)), active_rn))
    if(length(active_rn)==(nrow(reactions_full)-1)){
      reaction=setdiff((1:nrow(reactions_full)), active_rn)
    }
    else{
      reaction=sample(setdiff((1:nrow(reactions_full)), active_rn), 1)
    }
    active_rn=c(active_rn, reaction)
    #print(active_rn)
    
    reactions=reactions_full[active_rn,]
    #metabolites=performReactions(reactions[active_rn,], metabolites, SpeciesConst = SpeciesConst, stopCriteria = "Eqb")
    
    S=matrixFromList(reactions, metabolites)
    metabolites=performReactions_ode(metabolites, reactions, SpeciesConst = SpeciesConst, S = S)
    criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  }
  return(list(metabolites=metabolites, criticalPoints=criticalPoints, sequence=active_rn, profile=df))
}
'
thaw_oneByone_ode<-function(metabolites, reactions, SpeciesConst=NULL, print = F, detailed = F, maxConc=2){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(reactions, metabolites)
  }
  active_rn=NULL
  criticalPoints=NULL
  criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  reactions_full=reactions
  while(length(active_rn)<nrow(reactions_full)){
    #switching on rn randomly one by one
    #print(setdiff((1:nrow(reactions_full)), active_rn))
    if(length(active_rn)==(nrow(reactions_full)-1)){
      reaction=setdiff((1:nrow(reactions_full)), active_rn)
    }
    else{
      reaction=sample(setdiff((1:nrow(reactions_full)), active_rn), 1)
    }
    active_rn=c(active_rn, reaction)
    #print(active_rn)
    
    reactions=reactions_full[active_rn,]
    #metabolites=performReactions(reactions[active_rn,], metabolites, SpeciesConst = SpeciesConst, stopCriteria = 'Eqb')
    
    S=matrixFromList(reactions, metabolites)
    metabolites=performReactions_ode(metabolites, reactions, SpeciesConst = SpeciesConst, S = S,detailed=detailed, maxConc=maxConc)
    criticalPoints=c(criticalPoints, as.numeric(colnames(metabolites)[ncol(metabolites)]))
  }
  return(list(metabolites=metabolites, criticalPoints=criticalPoints, sequence=active_rn))
}

'
performReactions_ode<-function(metabolites, reactions, SpeciesConst=NULL, S=NULL, display=F){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(reactions=reactions, metabolites=metabolites)
  }
  if(is.null(S)){
    S=matrixFromList(reactions=reactions, metabolites=metabolites)
  }
  yini<-metabolites[,ncol(metabolites)]
  out1<-runsteady(time=c(0,1e5), y=yini, func=fluxCalc, 
                  parms=list(S, SpeciesConst, maxConc, reactions$rateConst))
  time=attr(out1, "time")
  dt=0.01 
  if(time/dt>10000){
    dt=time/10000
  }
  time_init=0
  if(ncol(metabolites)>3){
    time_init=as.numeric(colnames(metabolites)[ncol(metabolites)])
  }
  if(dt==0.01){
    times<-round(seq(from=time_init, to=time_init+attr(out1, "time")+dt*1000, by=dt),2)
  }
  else{
    times<-round(seq(from=time_init, to=time_init+attr(out1, "time"), by=dt),2)
  }
  out<-ode(times=times, y=yini, func=fluxCalc, 
           parms=list(S, SpeciesConst, maxConc, reactions$rateConst))
  if(ncol(metabolites)>3){
    out<-out[-1,]
  }
  colnames<-c(colnames(metabolites), as.character(out[,1]))
  out<-out[,-1]
  metabolites<-cbind(metabolites, t(out))
  colnames(metabolites)<-colnames
  
  if(display==T){
    plotConcProfile(metabolites = metabolites, log="x")
  }
  return(metabolites)
}
'

performReactions_ode<-function(metabolites, reactions, SpeciesConst=NULL, S=NULL, display=F, detailed=F, maxConc=2){
  if(is.null(SpeciesConst)){
    SpeciesConst=findConstSpecies(reactions=reactions, metabolites=metabolites)
  }
  if(is.null(S)){
    S=matrixFromList(reactions=reactions, metabolites=metabolites)
  }
  yini<-metabolites[,ncol(metabolites)]
  out1<-runsteady(time=c(0,1e5), y=yini, func=fluxCalc, 
                  parms=list(S, SpeciesConst, maxConc, reactions$rateConst))
  time=attr(out1, "time")
  if(detailed == F){
    colnames<-c(colnames(metabolites), as.character(round(time,2)))
    metabolites<-cbind(metabolites, out1$y)
    colnames(metabolites)<-colnames
  }else{
    dt=0.01 
    if(time/dt>10000){
      dt=time/10000
    }
    time_init=0
    if(ncol(metabolites)>3){
      time_init=as.numeric(colnames(metabolites)[ncol(metabolites)])
    }
    times<-round(seq(from=time_init, to=time_init+attr(out1, "time"), by=dt),2)
    out<-ode(times=times, y=yini, func=fluxCalc, 
             parms=list(S, SpeciesConst, maxConc, reactions$rateConst))
    if(ncol(metabolites)>3){
      out<-out[-1,]
    }
    colnames<-c(colnames(metabolites), as.character(out[,1]))
    out<-out[,-1]
    metabolites<-cbind(metabolites, t(out))
    colnames(metabolites)<-colnames
  }
  if(display==T){
    plotConcProfile(metabolites = metabolites, log="x", maxConc = maxConc)
  }
  return(metabolites)
}

source('data_collector.R')


metRxns_from_stoichMatrix<-function(stoich_matrix, lowerLimitConc=0, upperLimitConc=1){
  metabolites=data.frame(matrix(0,ncol(stoich_matrix),3))
  colnames(metabolites)=c('name', 'mass', 'conc')
  metabolites$name=1:ncol(stoich_matrix)
  metabolites$conc=round((runif(length(metabolites$name), min = lowerLimitConc, max = upperLimitConc))*10000)/10000
  
  #make the reaction matrix
  colnames=c('reactant1', 'reactant2', 'product1', 'product2', 'rateConst', 'concReactant1', 'concReactant2', 'rate')
  reactions=data.frame(matrix(0, nrow(stoich_matrix), 8))
  colnames(reactions)=colnames
  for(i in 1:nrow(stoich_matrix)){
    #print(i)
    reaction=stoich_matrix[i,]
    reactants=which(reaction<0)
    products=which(reaction>0)
    #print(reactants)
    #print(products)
    reactions[i, 1:length(reactants)]=reactants
    reactions[i, 3:(2+length(products))]=products
    reactions$rateConst=1
  }
  return(list(metabolites=metabolites, reactions=reactions))
}


findExtSpecies<-function(metabolites, reactions, metaboliteNet=NULL, lowerLimitConc=0, upperLimitConc=1){
  #In this function, I find which metabolites are not products in a reaction and which metabolites are not reactants in any reaction.
  #Such metabolites must be imported to and exported from the cell respectively.
  #For each metabolite, we add one metabolite which is the extracellular form of that metabolite. 
  #We will return the new metabolites and reactions Dataframe along with a list of the names of the external metabolites
  
  if(is.null(metaboliteNet)){
    metaboliteNet=getNetwork(metabolites = metabolites, reactions = reactions, print = F, display = F)$metaboliteNet
  }
  dg=decompose.graph(metaboliteNet)
  SpeciesConst=NULL
  sources=NULL
  sinks=NULL
  sources_tot=NULL
  sinks_tot=NULL
  for(i in 1:length(dg)){
    component=dg[[i]]
    if(length(vertex_attr(component)$name)==1){
      SpeciesConst=c(SpeciesConst, vertex_attr(component)$name)
      next
    }
    #find the indegree of each node of component
    #find the outdegree of each node of component
    #if none is zero for indegree then randoml choose some metbaoite
    #if none is zero for outdegree then reandomly choose some metabolite which is not the earlier one
    indegrees=degree(component, v=V(component), mode='in')
    sources=vertex_attr(component)$name[indegrees==0]
    outdegrees=degree(component, v=V(component), mode='out')
    sinks=vertex_attr(component)$name[outdegrees==0]
    if(length(sources)==0){
      #randomly choose one metabolite which is not in sinks
      sources=sample(setdiff(vertex_attr(component)$name, sinks), 1) 
    }
    if(length(sinks)==0){
      #randomly choose one metabolite which is not in sources
      sinks=sample(setdiff(vertex_attr(component)$name, sources), 1)
    }
    sources_tot=union(sources, sources_tot)
    sinks_tot=union(sinks, sinks_tot)
    SpeciesConst=c(SpeciesConst, sources, sinks)
  }
  SpeciesConst=sort(as.numeric(SpeciesConst))
  sources=sort(as.numeric(sources_tot))
  sinks=sort(as.numeric(sinks_tot))
  
  #metabolites which are neither sources nor sinks are isolated metabolites
  #Not adding external metabolites for isolated metabolites
  isolated_metabolites=setdiff(SpeciesConst, union(sources, sinks))
  SpeciesConst=intersect(SpeciesConst, union(sources, sinks))
  
  lastIndex_metabolites=nrow(metabolites)
  lastIndex_reactions=nrow(reactions)
  #Add the corresponding metabolites
  NbSpeciesConst=length(SpeciesConst)
  metabolites[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)),]=0
  externalMetabolites=(lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)
  metabolites$name[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=(lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)
  metabolites$mass[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=metabolites$mass[SpeciesConst]
  metabolites$conc[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=round((runif(length(SpeciesConst), min = lowerLimitConc, max = upperLimitConc))*10000)/10000
  
  
  reactions[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst), ]=0
  reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=SpeciesConst
  reactions$product1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=externalMetabolites
  #those reactions in which reactant1 is a source, their respective external metabolites are the reactants
  indices=which(reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)] %in% sources)+lastIndex_reactions
  reactions[indices, c(1,3)]=reactions[indices, c(3,1)]
  #updating the concentrations of the reactants.
  reactions$concReactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=metabolites$conc[reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]]
  reactions$concReactant2[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=rep(1, NbSpeciesConst)
  reactions$rateConst[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=rep(1, NbSpeciesConst)
  #Add the corresponding reactions
  
  return(list(metabolites=metabolites, reactions=reactions, externalMetabolites=union(externalMetabolites, isolated_metabolites)))
}


add_sources<-function(metabolites, reactions, SpeciesConst, sourceOf, lowerLimitConc=0, upperLimitConc=1, g=NULL){
  #This function adds new source metabolites for 'sourceOf' metabolites
  #This will not change the initial conc of the metabolites which were already existing 
  #this will also add reactions where new metabolites convert to 'sourceOf' metabolites
  #this will return new metabolites, reactions, and new external Species
  
  lastIndex_metabolites=nrow(metabolites)
  lastIndex_reactions=nrow(reactions)
  
  if(is.null(g)){
    g=getNetwork(metabolites = metabolites, reactions = reactions, SpeciesConst = SpeciesConst)
  }
  #find if any of the sourceOf metabolites already has a source metabolite
  #get the speciesConst of the network
  #check if the sourceOf metabolites are the neighbors of these external metabolites
  neighbors_ext=as.numeric(unlist(lapply(SpeciesConst, function(x) as_ids(neighbors(g$metaboliteNet, as.character(x))))))#neighbors of external metabolites
  if(length(which(sourceOf %in% neighbors_ext))>0){
    sourceOf=sourceOf[-which(sourceOf %in% neighbors_ext)]
  }
  
  
  NbSpeciesConst=length(sourceOf) # I am using the variable named NbSpeciesConst because 
  #this variable was being used in another code from which the next snippet has been taken.
  #Dont confuse this with length of SpeciesConst
  metabolites[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)),]=0
  externalMetabolites=(lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)
  
  metabolites$name[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=(lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)
  metabolites$mass[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=metabolites$mass[sourceOf]
  metabolites$conc[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=round((runif(NbSpeciesConst, min = lowerLimitConc, max = upperLimitConc))*10000)/10000
  
  reactions[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst), ]=0
  reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=externalMetabolites
  reactions$product1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=sourceOf
  reactions$concReactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=metabolites$conc[reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]]
  reactions$concReactant2[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=rep(1, NbSpeciesConst)
  reactions$rateConst[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=rep(1, NbSpeciesConst)
  
  externalSpecies=union(SpeciesConst, externalMetabolites)
  
  return(list(metabolites=metabolites, reactions=reactions, SpeciesConst=externalSpecies))
}

add_sinks<-function(metabolites, reactions, SpeciesConst, sinkOf, lowerLimitConc=0, upperLimitConc=1, g=NULL){
  #This function adds new source metabolites for 'sinkOf' metabolites
  #This will not change the initial conc of the metabolites which were already existing 
  #this will also add reactions where new metabolites convert to 'sinkOf' metabolites
  #this will return new metabolites, reactions, and new external Species
  
  lastIndex_metabolites=nrow(metabolites)
  lastIndex_reactions=nrow(reactions)
  
  if(is.null(g)){
    g=getNetwork(metabolites = metabolites, reactions = reactions, SpeciesConst = SpeciesConst)
  }
  #find if any of the sinkOf metabolites already has a source metabolite
  #get the speciesConst of the network
  #check if the sinkOf metabolites are the neighbors of these external metabolites
  neighbors_ext=as.numeric(unlist(lapply(SpeciesConst, function(x) as_ids(neighbors(g$metaboliteNet, as.character(x), mode='in')))))#neighbors of external metabolites
  if(length(which(sinkOf %in% neighbors_ext))>0){
    sinkOf=sinkOf[-which(sinkOf %in% neighbors_ext)]
  }
  
  NbSpeciesConst=length(sinkOf) # I am using the variable named NbSpeciesConst because 
  #this variable was being used in another code from which the next snippet has been taken.
  #Dont confuse this with length of SpeciesConst
  metabolites[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)),]=0
  externalMetabolites=(lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)
  
  metabolites$name[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=(lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst)
  metabolites$mass[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=metabolites$mass[sinkOf]
  metabolites$conc[((lastIndex_metabolites+1):(lastIndex_metabolites+NbSpeciesConst))]=round((runif(NbSpeciesConst, min = lowerLimitConc, max = upperLimitConc))*10000)/10000
  
  reactions[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst), ]=0
  reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=sinkOf
  reactions$product1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=externalMetabolites
  reactions$concReactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=metabolites$conc[reactions$reactant1[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]]
  reactions$concReactant2[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=rep(1, NbSpeciesConst)
  reactions$rateConst[(lastIndex_reactions+1):(lastIndex_reactions+NbSpeciesConst)]=rep(1, NbSpeciesConst)
  
  externalSpecies=union(SpeciesConst, externalMetabolites)
  
  return(list(metabolites=metabolites, reactions=reactions, SpeciesConst=externalSpecies))
}

