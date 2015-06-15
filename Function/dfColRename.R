#function that adds a part to the colum names of a data frame for identification,txtn = identifier POS=position to start
dfColRename<-function(df,txtn,POS){
  for (i in 1:ncol(df)){
    names(df)[i]=paste(txtn,substring(names(df)[i],POS),sep="")
  }
  return(df)
}
