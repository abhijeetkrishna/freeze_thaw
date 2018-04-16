#analyse a real network from kegg
#this code uses findExtSpecies instead of findConstSpecies

#pathway=
detailed = F
upperLimitConc<-1
lowerLimitConc<-0
maxConc=5
Nbrandom=10

source('functions.R')

set.seed=1
dir=paste0('kegg/data_files/')
stoich_matrix=read.csv(paste0(dir, 'stoich_matrix.csv'), header = F, stringsAsFactors = F)

name <- basename(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])
dir.create(paste0('kegg/',substr(name, 1, nchar(name)-2)))
pathway_id=00500
resultDir=paste0(paste0('kegg/',substr(name, 1, nchar(name)-2), '/', as.character(pathway_id),'/'))
dir.create(resultDir)

NbSpecies=ncol(stoich_matrix)
NbReactions=nrow(stoich_matrix)

print(NbSpecies)
print(NbReactions)

colnames<-c('case', 'i', as.character(1:NbSpecies)) #because of this we do not see the external species concentration
df=data.frame(matrix(0, 0, length(colnames)))
colnames(df)=colnames


seeds=1:10
for(seed in seeds){
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

