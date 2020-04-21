args <- commandArgs()
species <- args[5]

### 1.KO文件下载
## 物种转换，后续可自定义添加
time <- gsub("-","",Sys.Date())
species_ko <- tolower(species)
species_list <- data.frame(matrix(c("mouse","mmu","human","hsa"),ncol=2,byrow=T),stringsAsFactors=F)
names(species_list) <- c("Species","Species_kegg")
if(species_ko%in%species_list[,1]){
		species_ko <- species_list[grep(species_ko,species_list[,1]),2]
}


if(!length(grep(paste(species_ko,"_kegg.ko",sep=""),list.files()))>0){
		download.file(paste("https://www.genome.jp/kegg-bin/download_htext?htext=",species_ko,"00001&format=htext&filedir=",sep=""),destfile=paste(species_ko,"_kegg.ko_",time,sep=""))
		kegg_ko <- read.table(paste(species_ko,"_kegg.ko_",time,sep=""),header=F,sep="\t",fill=T,quote="",stringsAsFactors=F)
} else {
		kegg_ko <- read.table(list.files()[grep(paste(species_ko,"_kegg.ko_",sep=""),list.files())],header=F,sep="\t",fill=T,quote="",stringsAsFactors=F)
}
