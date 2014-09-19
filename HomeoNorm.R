#Simon Renny-Byfield, Iowa State University September 2014

#A function to normalize libraries with homeologous counts.
#The script assumes that each individual library has been 
#passed through polyCat and been split into At,Dt and N counts.
#Each sequence of three columns will be treated as a single library
#and the normalization number for each will be calculated based 
#on the total counts (At, Dt and N) from the library (now in three seperate 
#columns). This approach ensures the ratios betwee At, Dt and N are maintained
#within a sample, even after normalization between samples has
#occured.

#######
#A function to return a normalized read count table, in RPM
#######

HomeoNorm<- function ( count ) {
  #save the gene names
  names<-count[,1]
  #the dimensions of the table assuming the first column is gene names
  count<-count[,-1]
  dims<-dim(count)
  #set up a new data.frame to contain the modifed counts
  NormCounts<-count
  #check to see if the number of columns is a multiple of 3
  if ( dims[2] %% 3 == 0 ) {
    #find the number of libraries
    NumLib<-dims[2]/3
    #a sequence of starting columns
    ColSeq<-seq(from=1, to=dims[2], by =3)
    #print(ColSeq)
    libSizes<-list()
    #figure out each library size
    for ( i in 1:length(ColSeq) ) {
      print(i)
      #the sequence of columns to grab
      toModify<-seq(from=ColSeq[i], to=ColSeq[i]+2)
      lib<-count[,toModify]
      #calculate the tota library size
      libSizes[[i]]<-sum(colSums(lib))
      #normalize all three columns, At, Dt and N, based on the total lib size  
      NormCounts[,toModify]<-sweep(NormCounts[,toModify],2,libSizes[[i]],FUN="/")
      NormCounts[,toModify]<-NormCounts[,toModify]*1e6
    }#for    
  }#if
  else {
    stop('Incorrect number of columns (not a multiple of three) \n\n')
  }#else
  #return the norm counts and library sizes
  output<-list("normalized"=NormCounts,"sizes"=libSizes)
  return(output)
}#NormFunction
