
chain_metabolites<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  graph=g[[3]]
  chains=NULL
  lengths=NULL
  vertices=vertex_attr(graph)$name
  vertices1=vertices[unlist(lapply(vertices, function(x) 
    (degree(graph, x, mode = 'in')<=1 & degree(graph, x, mode='out')<=1)))]
  linearGraph=induced_subgraph(graph, vertices1)
  components=decompose(linearGraph, mode = "weak", max.comps = NA,
                       min.vertices = 0)
  if(length(components)==0){
    return(list(maxlength=0, sumlength=0, Nbchains=0))
  }
  for(i in 1:length(components)){
    #find the component
    component=components[[i]]
    #if the size of component <=1 next
    if(length(as_ids(V(component)))<=1) next
    #Make sure that the component is not a loop
    if(min(degree(component, V(component)), mode='all')==2) next
    #Add it to chains
    chains[[length(chains)+1]]=as_ids(V(component))
    #update lengths 
    lengths=c(lengths, length(as_ids(V(component))))
  }
  return(list(maxlength=ifelse(max(lengths)==-Inf, 0, max(lengths)), sumlength=sum(lengths), Nbchains=length(lengths)))
}

chain_reactions<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  graph=g[[2]]
  chains=NULL
  lengths=NULL
  vertices=vertex_attr(graph)$name
  vertices1=vertices[unlist(lapply(vertices, function(x) 
    (degree(graph, x, mode = 'in')<=1 & degree(graph, x, mode='out')<=1)))]
  linearGraph=induced_subgraph(graph, vertices1)
  components=decompose(linearGraph, mode = "weak", max.comps = NA,
                       min.vertices = 0)
  if(length(components)==0){
    return(list(maxlength=0, sumlength=0, Nbchains=0))
  }
  for(i in 1:length(components)){
    #find the component
    component=components[[i]]
    #if the size of component <=1 next
    if(length(as_ids(V(component)))<=1) next
    #Make sure that the component is not a loop
    if(min(degree(component, V(component)), mode='all')==2) next
    #Add it to chains
    chains[[length(chains)+1]]=as_ids(V(component))
    #update lengths 
    lengths=c(lengths, length(as_ids(V(component))))
  }
  return(list(maxlength=ifelse(max(lengths)==-Inf, 0, max(lengths)), sumlength=sum(lengths), Nbchains=length(lengths)))
}

Nb_branching_metabolites<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #count the number of metabolites having outdegree > 1
  graph=g[[3]]
  degrees=degree(graph, V(graph), mode = 'out')
  vertices=as_ids(V(graph))
  return(length(vertices[degrees>=2]))
}

Nb_branching_reactions<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #count the number of reactions having outdegree > 1
  graph=g[[3]]
  degrees=degree(graph, V(graph), mode = 'out')
  vertices=as_ids(V(graph))
  return(length(vertices[degrees>=2]))
}

Nb_feedback_reactions<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #for each pair of nodes, check if there exists a simple path from 1 to 2
  #If the path exists then check if there exists path from 2 to 1
  graph=g[[2]]
  vertices=as_ids(V(graph))
  list=NULL
  lengths=NULL
  for(vertex1 in vertices){
    for(vertex2 in vertices){
      simple_paths=all_simple_paths(graph, from = vertex1, to = vertex2, mode='out')
      if(length(simple_paths)==0) next
      if(get.edge.ids(graph, vp=c(vertex2,vertex1), directed = T)>0){
        #add this list to the list of loops
        list=c(list, simple_paths)
        lengths=c(lengths, unlist(lapply(1:length(simple_paths), function(x) length(simple_paths[[x]]))))
      }
    }
  }
  #lengths=c(2,3,4,3,2,3,3,4,3,3,2,2,4,2,2,4)
  
  lengths=lengths[lengths>2]
  
  size_loops=0
  freq_loops=0
  if(length(lengths)>0){
    l=split(lengths, lengths)
    size_loops=unlist(lapply(1:length(l), function(x) unique(l[[x]])))
    freq_loops=unlist(lapply(1:length(l), function(x) length(l[[x]])/unique(l[[x]])))
  }
  
  return(list(nb_loops=sum(freq_loops), size_loops=size_loops, freq_loops=freq_loops))
}

Nb_feedforward_reactions<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  graph=g[[2]]
  vertices=as_ids(V(graph))
  list=NULL
  lengths=NULL
  for(vertex1 in vertices){
    for(vertex2 in vertices){
      simple_paths=all_simple_paths(graph, from = vertex1, to = vertex2, mode='out')
      if(length(simple_paths)<=1) next
      if(get.edge.ids(graph, vp=c(vertex1,vertex2), directed = T)>0){
        #add this list to the list of loops
        list=c(list, simple_paths)
        lengths=c(lengths, unlist(lapply(1:length(simple_paths), function(x) length(simple_paths[[x]]))))
      }
    }
  }
  #lengths=c(2,3,4,3,2,3,3,4,3,3,2,2,4,2,2,4)
  lengths=lengths[lengths>2]
  size_loops=0
  freq_loops=0
  if(length(lengths)>0){
    size_loops=sort(unique(lengths))
    l=split(lengths, lengths)
    freq_loops=unlist(lapply(1:length(l), function(x) length(l[[x]])))
  }
  return(list(nb_loops=sum(freq_loops), size_loops=size_loops, freq_loops=freq_loops))
}

Nb_feedback_metabolites<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  graph=g[[3]]
  vertices=as_ids(V(graph))
  list=NULL
  lengths=NULL
  for(vertex1 in vertices){
    for(vertex2 in vertices){
      simple_paths=all_simple_paths(graph, from = vertex1, to = vertex2, mode='out')
      if(length(simple_paths)==0) next
      if(get.edge.ids(graph, vp=c(vertex2,vertex1), directed = T)>0){
        #add this list to the list of loops
        list=c(list, simple_paths)
        lengths=c(lengths, unlist(lapply(1:length(simple_paths), function(x) length(simple_paths[[x]]))))
      }
    }
  }
  #lengths=c(2,3,4,3,2,3,3,4,3,3,2,2,4,2,2,4)
  lengths=lengths[lengths>2]
  size_loops=0
  freq_loops=0
  if(length(lengths)>0){
    l=split(lengths, lengths)
    size_loops=unlist(lapply(1:length(l), function(x) unique(l[[x]])))
    freq_loops=unlist(lapply(1:length(l), function(x) length(l[[x]])/unique(l[[x]])))
  }
  
  return(list(nb_loops=sum(freq_loops), size_loops=size_loops, freq_loops=freq_loops))
}

Nb_feedforward_metabolites<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  graph=g[[3]]
  vertices=as_ids(V(graph))
  list=NULL
  lengths=NULL
  for(vertex1 in vertices){
    for(vertex2 in vertices){
      simple_paths=all_simple_paths(graph, from = vertex1, to = vertex2, mode='out')
      if(length(simple_paths)<=1) next
      if(get.edge.ids(graph, vp=c(vertex1,vertex2), directed = T)>0){
        #add this list to the list of loops
        list=c(list, simple_paths)
        lengths=c(lengths, unlist(lapply(1:length(simple_paths), function(x) length(simple_paths[[x]]))))
      }
    }
  }
  #lengths=c(2,3,4,3,2,3,3,4,3,3,2,2,4,2,2,4)
  lengths=lengths[lengths>2]
  
  size_loops=0
  freq_loops=0
  if(length(lengths)>0){
    size_loops=sort(unique(lengths))
    l=split(lengths, lengths)
    freq_loops=unlist(lapply(1:length(l), function(x) length(l[[x]])))
  }
  
  return(list(nb_loops=sum(freq_loops), size_loops=size_loops, freq_loops=freq_loops))
}

Nb_met_NCTEM<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #metabolites not connected to external metabolites
  graph=g[[3]]
  vertices=as_ids(V(graph))
  #label each vertex as connected or not-connected
  connected=NULL
  notConnected=NULL
  for(vertex in vertices){
    if(vertex %in% as.character(SpeciesConst)){
      connected=c(connected, vertex)
      next
    } 
    flag=0
    for(extMet in SpeciesConst){
      paths=all_simple_paths(graph, from = vertex, to = as.character(extMet))
      if(length(paths)>0){
        flag=1
        break
      }
    }
    if(flag==1){
      connected=c(connected, vertex)
    }else{
      notConnected=c(notConnected, vertex)
    }
  }
  return(list(NbNotConnected=length(notConnected), connected=connected, notConnected=notConnected))
}

Nb_reversibleRxn<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #for each reaction find if the opposite exists
  #might be faster throgh stoichiometric matrix
  #Make a new matrix where each rxn is repeated for nrow(reactions) times
  #Make a matrix where reactions is repeated for 
  graph=g[[2]]
  S=matrixFromList(reactions, metabolites)
  S=t(S) 
  ####### Note that I am using S in a way that each row is reaction and each col is metabolite
  avoid=NULL
  df_rev=data.frame(matrix(0, 0, 2))
  counter=0
  for(i in 1:nrow(reactions)){
    if(i %in% avoid) next
    rxn=S[i,]
    rxnMat=matrix(rxn, nrow(S), ncol(S), byrow=T)
    M=rxnMat+S
    #find index of row which is zero
    index=which(unlist(lapply(1:nrow(reactions), function(x) all(M[x,]==rep(0,nrow(metabolites))))))
    if(length(index)>0){
      avoid=c(avoid, index)
      counter=counter+1
      df_rev[counter,]=c(i, index)
    } 
  }
  return(nrow(df_rev))
}

'Nb_revIrrevChain<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
}'

Nb_magnifyGlass_reactions<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #Find all the single loop reactions
  #Check the outdegrees and indegrees of the two reactions
  #Subtract 2 from total outdegrees and indegrees
  #One of these quantities should be zero and the other should be non-zero
  #also returns the number of simple loops
  graph=g[[2]]
  count=0
  Nbsimpleloops=0
  #i=0
  edgelist=as_edgelist(graph)
  while(nrow(edgelist)>1){
    #if(i==12) break
    #print(nrow(edgelist))
    #i=i+1
    #print(i)
    index=which(unlist(lapply(2:nrow(edgelist), function(x) all(edgelist[x,]==edgelist[1,c(2,1)]))))+1
    if(length(index)==1){
      Nbsimpleloops=Nbsimpleloops+1
      vertices=edgelist[1,]
      outdegree=sum(degree(graph, vertices, mode='out'))
      indegree=sum(degree(graph, vertices, mode='in'))
      outdegree=outdegree-2
      indegree=indegree-2
      if(outdegree==0 | indegree==0){
        count=count+1
      }
      edgelist=edgelist[-index,]
      
    }
    if(is.null(nrow(edgelist))) break #this is used when edgelist has only one row
    if(nrow(edgelist)<=2) break
    edgelist=edgelist[-1,]
  }
  return(list(countMagnifyGlass=count, Nbsimpleloops=Nbsimpleloops))
}

Nb_magnifyGlass_metabolites<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  #Find all the single loop reactions
  #Check the outdegrees and indegrees of the two reactions
  #Subtract 2 from total outdegrees and indegrees
  #One of these quantities should be zero and the other should be non-zero
  #also returns the number of simple loops
  graph=g[[3]]
  count=0
  Nbsimpleloops=0
  #i=0
  edgelist=as_edgelist(graph)
  while(nrow(edgelist)>1){
    #if(i==12) break
    #print(nrow(edgelist))
    #i=i+1
    #print(i)
    index=which(unlist(lapply(2:nrow(edgelist), function(x) all(edgelist[x,]==edgelist[1,c(2,1)]))))+1
    if(length(index)==1){
      Nbsimpleloops=Nbsimpleloops+1
      vertices=edgelist[1,]
      outdegree=sum(degree(graph, vertices, mode='out'))
      indegree=sum(degree(graph, vertices, mode='in'))
      outdegree=outdegree-2
      indegree=indegree-2
      if(outdegree==0 | indegree==0){
        count=count+1
      }
      edgelist=edgelist[-index,]
      
    }
    if(is.null(nrow(edgelist))) break #this is used when edgelist has only one row
    if(nrow(edgelist)<=2) break
    edgelist=edgelist[-1,]
  }
  return(list(countMagnifyGlass=count, Nbsimpleloops=Nbsimpleloops))
}

'Nb_genMagnifyGlass<-function(g, metabolites=NULL, reactions=NULL, SpeciesConst=NULL){
  
}'


# From all the loops, remove those loops which have species const in their loops



