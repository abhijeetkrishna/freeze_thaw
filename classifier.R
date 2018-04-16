
source('functions.R')
dir.create('classified')
path= 'Thaw_oneByone_extMet/'
folders=list.files(path)
Nbi=10
for(folder in folders){
  NbSpecies=as.numeric(substr(folder, 1, regexpr('metabolites', folder)-1))
  NbReactions=as.numeric(substr(folder, regexpr('_', folder)+1, regexpr('reactions', folder)-1))
  #NbSpecies=4 #if you want to classify a particular folder
  #NbReactions=3 #if you want to classify a particular folder
  path= 'Thaw_oneByone_extMet/'
  path= paste0(path, as.character(NbSpecies), 'metabolites_', as.character(NbReactions), 'reactions/')
  resultDir = paste0('classified/',as.character(NbSpecies), 'metabolites_', as.character(NbReactions), 'reactions/')
  dir.create(resultDir)
  dir.create(paste0(resultDir,'Sequence_Dependent/'))
  dir.create(paste0(resultDir,'Sequence_Independent/'))
  dir.create(paste0(resultDir,'Thaw_Independent/'))
  files_forSeed=list.files(path, 'thawedProfile_seed_')
  count=0
  for(file in files_forSeed){
    file=paste0(path, file)
    count=count+1
    print( paste0('NbSpecies = ', as.character(NbSpecies), ' NbReactions = ', as.character(NbReactions), ' count = ', as.character(count)))
    seed=as.numeric(substr(file, regexpr('seed_', file)+nchar('seed_'), regexpr('.csv', file)-1))
    if(!file.exists(file)) next
    df=read.csv(file, stringsAsFactors = F)
    mat1=as.numeric(matrix(df[1,4:ncol(df)], nrow(df)-1, ncol(df)-3, byrow = T))
    mat2=as.numeric(matrix(df[2,4:ncol(df)], nrow(df)-1, ncol(df)-3, byrow = T))
    mat3=as.numeric(as.matrix(df[2:nrow(df), 4:ncol(df)]))
    mat4=mat1-mat3
    mat5=mat2-mat3
    files=NULL
    files=c(file, paste0(path,'network_seed_', as.character(seed), '.pdf'))
    files=c(files, paste0(path, 'profile_seed_', as.character(seed), '_i_',
                          as.character(1:Nbi), '.csv'))
    if(any(abs(mat5)>0.1)){
      #copy the files to the folder sequence dependent
      resultDir1<-paste0(resultDir,'Sequence_Dependent/')
      
    }else if(any(abs(mat4)>0.1)){
      #copy the files to the folder sequence independent
      resultDir1<-paste0(resultDir,'Sequence_Independent/')
    }else{
      #copy the files to the folder single state
      resultDir1<-paste0(resultDir,'Thaw_Independent/')
    }
    for(j in 1:length(files)){
      file.copy(files[j], resultDir1)
    }
  }
  #break #if you want to classify a particular folder
}
