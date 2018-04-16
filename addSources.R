#This code identifies the metabolites which are most perturbed and adds transporters to them 


#Identifying the metabolites which are most perturbed. 
source('functions.R')
pathway_id=00500
path=paste0('kegg/kegg_analysis_extMet/', as.character(pathway_id),'/')
resultDir=path
seed=1
maxConc=5
lowerLimitConc=0
upperLimitConc=1
detailed=F
Nbrandom=10

files=list.files(path, 'thawedProfile_')

flag=0

for(file in files){
  #file=files[1]
  print(file)
  temp=read.csv(paste0(path, file))
  temp=temp[,4:ncol(temp)]
  
  if(flag==0){
    flag=1
    df=data.frame(matrix(0, 0, ncol(temp)))
    #once you get the name of all the metabolites, make the names the colnames
  }
  
  eqb=temp[1,]
  for(i in 2:nrow(temp)){#first one is just equilibrium concentraitons
    diff=temp[i,]-eqb
    df[nrow(df)+1,]=diff
    #for(j in 1:ncol(temp)){
    #subtract and add it in the correct place
    #}
  }
}

#write.csv(df, paste0(resultDir, 'diff_conc.csv'))
means=unlist(lapply(1:ncol(df), function(x) mean(df[,x])))
sds=unlist(lapply(1:ncol(df), function(x) sd(df[,x])))
lb_metabolites=means-2*sds
changing_metabolites=which((lb_metabolites)>0)


### Loading the real network
dir=paste0('kegg/data_files/')
stoich_matrix=read.csv(paste0(dir, 'stoich_matrix.csv'), header = F, stringsAsFactors = F)
NbReactions=nrow(stoich_matrix)
NbSpecies=ncol(stoich_matrix)

#Changing the resultDir to make sure that the data is stored at the correct place
resultDir='kegg/addSources/'
dir.create(resultDir)

seeds=1:10
for(seed in seeds){
  
  ################################################
  
  set.seed(seed)
  print(seed)
  list=metRxns_from_stoichMatrix(stoich_matrix)
  metabolites=list$metabolites
  reactions=list$reactions
  
  list2=findExtSpecies(metabolites, reactions)
  metabolites=list2$metabolites
  reactions=list2$reactions
  SpeciesConst=list2$externalMetabolites
  #SpeciesConst=findConstSpecies(metabolites, reactions)
  
  #Add external metabolites for the metabolites which are changing the most
  g=getNetwork(metabolites = metabolites, reactions = reactions, SpeciesConst = SpeciesConst)
  list3=add_sources(metabolites, reactions, SpeciesConst, sourceOf=changing_metabolites,g=g)
  metabolites=list3$metabolites
  reactions=list3$reactions
  SpeciesConst=list3$SpeciesConst
  
  ################################################
  
  #Make sure that while mapping the ids of the metabolites, the external metabolites do not raise an error
  
  #print(SpeciesConst)
  getNetwork(metabolites = metabolites, reactions = reactions,
             print = T, resultDir = resultDir, SpeciesConst = SpeciesConst,
             filename = paste0('network_seed_', as.character(seed), '.pdf'))
  #print(seed)
  S=matrixFromList(reactions, metabolites)
  metabolites=performReactions_ode(metabolites, reactions, SpeciesConst=SpeciesConst, S=S, detailed = detailed, maxConc = maxConc) #
  metabolites1=metabolites
  
  colnames<-c('case', 'i', as.character(1:NbSpecies))
  df=data.frame(matrix(0, 0, length(colnames)))
  colnames(df)=colnames
  
  df[1,]=c('RandomToEquilibrium', 0, metabolites[, ncol(metabolites)])
  for(i in 1:Nbrandom){
    #print(paste0('thaw one by one i = ', as.character(i)))
    #thawedProfile=thaw_oneByone_runsteady(metabolites=metabolites, eqbProfile=eqbProfile, reactions=reactions) #merge the two functions later
    print(paste0(as.character(i), 'th sequence'))
    metabolites<-metabolites1
    list2=thaw_oneByone_ode(metabolites=metabolites, reactions=reactions, SpeciesConst = SpeciesConst,detailed = detailed, maxConc=maxConc) #merge the two functions later
    metabolites<-list2[[1]]
    criticalPoints<-list2[[2]]
    sequence<-list2[[3]]
    
    df[i+1,]=c('Thawed', i, metabolites[, ncol(metabolites)])
    
    #outputs
    if(detailed==T){
      pdf(paste0(resultDir, 'thaw_oneByone_seed_', as.character(seed), '_i_', as.character(i), '.pdf'))
      plotConcProfile(metabolites, criticalPoints=criticalPoints, seed=seed, maxConc = maxConc)
      dev.off()
      write.csv(metabolites, paste0(resultDir, 'metabolites_seed_', as.character(seed), '_i_', as.character(i), '.csv'))
    }
    if(detailed==F){
      profile=metabolites
      colnames(profile)=c('name', 'mass', 'initialConc', 'Eqb', paste0('R', sequence))
      write.csv(profile, paste0(resultDir, 'profile_seed_', as.character(seed), '_i_', as.character(i), '.csv'))
    }
  }
  write.csv(df, paste0(resultDir, 'thawedProfile_seed_', as.character(seed), '.csv'))
  #seed=seed+1
}

