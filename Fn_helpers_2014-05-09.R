
# 
RangeFn = function(Vec) range(ifelse(abs(Vec)==Inf,NA,Vec),na.rm=TRUE) 

# transforms
f = function(Num){ if( var(as.vector(Num),na.rm=TRUE)==0 ){ Return = array(rep(0.5,prod(dim(Num))),dim(Num))}else{Return=(plogis(Num)-min(plogis(Num),na.rm=TRUE))/diff(range(plogis(Num),na.rm=TRUE))}; return(Return) } 
  #f = function(Num) array(order(Num),dim=dim(Num))/length(as.vector(Num))

# Colors
Col = colorRampPalette(colors=c("blue","purple","red"))

# Bins
Bin_Quantile = function(Obj, Nregions=4){
  Nbreaks = Nregions - 1
  Breaks = quantile(Obj, prob=seq(0,1,length=Nbreaks+2))
  Region = sapply( Obj, FUN=function(Num){ sum(Num>=Breaks) })
    Region = ifelse(Region==(Nregions+1),Nregions,Region)
    if( is.array(Obj) ) Region = array( Region, dim=dim(Obj))
  Return = list("Nregions"=Nregions, "Region"=Region, "Lwr"=Breaks[-c(Nbreaks+2)], "Upr"=Breaks[-1])
  return(Return)
}

# Make legend figure
Legend = function( Bin, Col, RowSet=NULL, Pcex=3, Tcex=3, Digits=3 ){
  if( is.null(RowSet) ) RowSet=1:Bin$Nregions
  plot( 1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame.plot=FALSE, xaxt="n", yaxt="n", xaxs="i", yaxs="i", mar=c(0,0,0,0) )
  Y = seq( 0, 1, length=length(RowSet)+2)[-c(length(RowSet)+1:2)] + 1/length(RowSet)/2
  points( x=rep(0.1,length(RowSet)), y=Y, col=Col(Bin$Nregions)[RowSet], pch=20, cex=Pcex)
  text( x=rep(0.2,length(RowSet)), y=Y, labels=paste(formatC(Bin$Lwr[RowSet],format="f",digits=Digits),"to",formatC(Bin$Upr[RowSet],format="f",digits=Digits)), pos=4, cex=Tcex)
}

# Bin = Bin_Quantile(Mat)
# Legend( Bin=Bin, Col=Col )
