#This is a thaw_oneByone_extMet version in which I am using findExtSpecies instead of findConstSpecies

source('functions.R')


##############
# Parameters #
##############

NbSpecies<-5
NbReactions<-5

minrange<-1
maxrange<-20
upperLimitConc<-1
lowerLimitConc<-0
maxConc=5
seed=1
Nbrandom=10
count=0
maxCount=400

detailed=F

name <- basename(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])
dir.create(substr(name, 1, nchar(name)-2))

for(k in 10:15){
  for(j in 10:15){
    
    NbSpecies=k
    NbReactions=j
    print(paste0('NbSpecies = ', as.character(NbSpecies), ' Nbreactions = ', as.character(NbReactions)))
    resultDir=paste0(substr(name, 1, nchar(name)-2), '/', as.character(NbSpecies), 
                     'metabolites_', as.character(NbReactions), 'reactions/')
    dir.create(resultDir)
    
    colnames<-c('case', 'i', as.character(1:NbSpecies))
    df=data.frame(matrix(0, 0, length(colnames)))
    colnames(df)=colnames
    #seed=1
    
    files=list.files(resultDir, pattern = 'network_seed_')
    if(length(files)>0){
      pos=unlist(lapply(1:length(files), function(x) regexpr('.pdf', files[x])[1]-1))
      seeds=as.numeric(unlist(lapply(1:length(files), function(x) substr(files[x], nchar('network_seed_')+1, pos[x]))))
      max_seed=max(seeds)
      seed=max_seed+1
    }else{
      seed=1
    }

    while(count<maxCount){
      set.seed(seed)
      list1<-initNetwork(NbSpecies, minrange, maxrange, upperLimitConc, lowerLimitConc, NbReactions=NbReactions)
      #list1<-initNetwork(NbSpecies, minrange, maxrange, upperLimitConc, lowerLimitConc, p=p)
      metabolites<-list1[[1]]
      reactions<-list1[[2]] 
      
      if(nrow(reactions)<NbReactions){
        seed=seed+1
        next
      }else{
        count=count+1
        print(paste0('count = ', as.character(count) , ' seed = ', as.character(seed)))
      }
      list2=findExtSpecies(metabolites, reactions)
      metabolites=list2$metabolites
      reactions=list2$reactions
      SpeciesConst=list2$externalMetabolites
      getNetwork(metabolites = metabolites, reactions = reactions,
                 print = T, resultDir = resultDir, SpeciesConst = SpeciesConst,
                 filename = paste0('network_seed_', as.character(seed), '.pdf'))
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
      seed=seed+1
    }
  }
}
