
args = commandArgs(T)

library(ggplot2)
library(reshape2)
library(Rmisc)
library(ggsci)

infile = args[1]
outfile = args[2]


#infile = "gene_exp20191219063Vi317.xls"
#outfile = "gene_exp20191219063Vi317.multi_plot.pdf"

con <- file(infile,'r')


comp_plot <- function(inline_data){
  
  header = inline_data[1]
  contont = inline_data[2]
  header_list = strsplit(header,"\t")[[1]]
  contont_list = strsplit(contont,"\t")[[1]]
#  print (contont_list)
  exp = rbind(header_list,contont_list)

  exp = as.data.frame(exp)
  colnames(exp) = header_list
  exp = exp[-1,]

  title0 <- strsplit(as.character(exp[1,1]),"/")[[1]][1]
  platform=strsplit(as.character(exp[1,1]),"/")[[1]][2]
  if( platform== "4.Comparison"){
    vs = sub(".total_genes.xls","",exp[1,1])[[1]]
    compare=strsplit(vs,"/")[[1]][3]
    sam1 = strsplit(compare,"_vs_")[[1]][1]
    sam2 = strsplit(compare,"_vs_")[[1]][2]
    
    grep1 <- grep(paste("^",sam1,".*(norm)",sep=""),colnames(exp))
    grep2 <- grep(paste("^",sam2,".*(norm)",sep=""),colnames(exp))
    
    
  }else if(platform == "2.DEG"){
    vs = sub("all_exp_","",exp[1,1])[[1]]
    vs = sub(".xls","",vs)
    compare=strsplit(vs,"/")[[1]][3]
    sam1 = strsplit(compare,"-vs-")[[1]][1]
    sam2 = strsplit(compare,"-vs-")[[1]][2]
    
    grep1 <- grep(paste("^",sam1,"_\\d+$",sep=""),colnames(exp))
    grep2 <- grep(paste("^",sam2,"_\\d+$",sep=""),colnames(exp))
    grep1_len=length(grep1)
    grep2_len=length(grep2)
    if((grep1_len ==0) |(grep2_len==0)){
        grep1 <- grep(paste("^",sam1,"_rep\\d+$",sep=""),colnames(exp))
        grep2 <- grep(paste("^",sam2,"_rep\\d+$",sep=""),colnames(exp))
    }
  }
  print(grep1)
  print(grep2)

  
  mydata <- melt(exp[,c(2,grep1,grep2)],id.vars=1)

  mydata$value=as.numeric(mydata$value)
  mean <- mean(as.numeric(as.character(unlist(exp[,grep1]))))+0.1
  mydata$value <- (mydata$value+0.1)/mean

  if(length(grep("TRUE",(mydata$value>2)))>0) {
    mydata[mydata$value>2,]$value <- 2 
  }
  #print(mydata)
  colnames(mydata) <- c("gene","sample","value")
  print(mydata$sample)
  if( platform== "4.Comparison"){
    mydata$stage <- sub("-[^-]+$","",mydata$sample)
  }else if(platform== "2.DEG"){
   contain_reps=grep("_rep\\d+$",mydata$sample)
   if(length(contain_reps)==0){
    mydata$stage <- sub("_\\d+$","",mydata$sample)
   }else{
    mydata$stage <- sub("_rep\\d+$","",mydata$sample)
   }
  }
  tgc <- summarySE(mydata,measurevar="value",groupvars=c("gene","stage"))
  tgc[is.na(tgc)] = 0.1
  tgc$stage <- factor(tgc$stage,levels=c(sam1,sam2))
  tgc$data <- exp[1,1]
  #tgc$stage <- as.character(tgc$stage)
  
  p = as.numeric(as.character(exp$pvalue))
  
  if(length(p)!=0) {
    if(p<0.05) {show_lable <- "*"; show_size <- 4}
    if(p<0.01) {show_lable <- "**"; show_size <- 4}
    if(p>0.05) {show_lable <- "n.s."; show_size <- 6/3}
  } else { show_lable <- "n.s."; show_size <- 6/3 }
  
  end1 <- tgc$value[1] + tgc$se[1]
  end2 <- tgc$value[2] + tgc$se[2]
  intercept1 <- (max(end1,end2))/15
  intercept2 <- (max(end1,end2))/10
  
  
  re_sam1=gsub("_","\n",sam1)
  re_sam2=gsub("_","\n",sam2)
  
  
  p1 <- ggplot(tgc,aes(x=stage,y=value,fill=stage)) +
    geom_bar(stat="identity",position="dodge",width=0.5,colour="black") +
    geom_errorbar(aes(ymin=value-se,ymax=value+se),colour="black",width=0.2) +
    geom_point(data=mydata,position=position_jitter(width=0.1,height=0),size=0.5,colour="black",aes(x=stage,y=value)) +
    geom_segment(aes(x=1,y=end1+intercept1,xend=1,yend=max(end1,end2)+intercept2)) +
    geom_segment(aes(x=2,y=end2+intercept1,xend=2,yend=max(end1,end2)+intercept2)) +
    geom_segment(aes(x=1,y=max(end1,end2)+intercept2,xend=1.3,yend=max(end1,end2)+intercept2)) + 
    geom_segment(aes(x=1.7,y=max(end1,end2)+intercept2,xend=2,yend=max(end1,end2)+intercept2)) +
    annotate("text",x=1.5,y=max(end1,end2)+intercept2,label=show_lable,size=show_size) +
    xlab(title0) +
    ylab("") +
    theme_bw() +
    theme(axis.text.x=element_text(size=6,color="black"),
          axis.title=element_text(size=6,color="black"),
          plot.title=element_text(size=6,color="black"),
          panel.grid=element_blank()) +
    theme(plot.margin=unit(c(0.05,0,0,0.05),"cm")) +
    guides(fill=FALSE,colour=FALSE)+
    ylim(0,2.5)+scale_x_discrete(labels=c(re_sam1,re_sam2))
  
  if(show_lable!="n.s.") {
    p2 <- p1 + theme(panel.backgr=element_rect(fill="#FFFF9B")) 
    return(p2)
  } else { 
  return (p1)
  }
}
i =1

while(TRUE){

  inline = readLines(con, n = 4)
  if(length(inline) == 0){
    break
  }
  com <- paste("plot",i," <- comp_plot(inline)",sep="")
  eval(parse(text=com))
  i=i+1
}

a=i-1


if (a<10){
h<-1.5
}else{
h<-1.3*ceiling(a/9)
}

fcom <- paste("multiplot(",paste("plot",1:a,sep="",collapse=","),",cols=10)",sep="")

pdf_plot <- paste("ggsave(",fcom,",filename='",outfile,"',width=13,height=",h,")",sep="")
eval(parse(text=pdf_plot)) 
