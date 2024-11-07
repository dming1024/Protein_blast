
#multiple blast of proteins
#包含有多条蛋白氨基酸序列的fasta文件
library(msa)
library(seqinr)
library(ape)

#这些基因的蛋白序列，是通过手动在uniproKB上获取的：https://www.uniprot.org/id-mapping

proteionBlast<-function(Gene){
  mySequenceFile <- sprintf("./fasta/%s.fa",Gene)
  mySequences <- readAAStringSet(mySequenceFile)
  #run the clustalW
  myFirstAlignment <- msa(mySequences)
  #输出pretty的结果
  msaPrettyPrint(myFirstAlignment, output="tex", file=sprintf("result_%s.tex",Gene),
                 showNames="left", showLogo="top",
                 shadingMode="identical",consensusColor="ColdHot",
                 #shadingModeArg = "standard area",
                 #furtherCode = c("\\showruler{1}{top}"),
                 showLegend=FALSE, askForOverwrite=FALSE)
  #设置miktex路径
  #Sys.getenv("PATH")
  #Sys.getenv('R_HOME')
  #Sys.which("pdflatex")
  #Sys.setenv(PATH=paste(Sys.getenv("PATH"),"D:\\proteinBlast0223\\miktex\\soft\\miktex\\bin\\x64",sep=";"))
  #修改tex文件中的texshade路径
  #"C:/Users/fan_qiangqiang/AppData/Local/Temp/"
  f=readLines(sprintf("result_%s.tex",Gene))
  f[[25]]=gsub("FAN_QI~1","fan_qiangqiang",f[[25]])
  f[[18]]=sprintf("\\topmargin=-%sin",1.5)
  writeLines(f,sprintf("result_%s.tex",Gene))
  tools::texi2pdf(sprintf("result_%s.tex",Gene),clean = T)
  
  
  hemoAln2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
  d <- dist.alignment(hemoAln2, "identity")
  hemoTree <- nj(d)
  pdf(sprintf("phylogenetic_of_%s.pdf",Gene))
  plot(hemoTree, main=sprintf("Phylogenetic Tree of %s Genes",Gene))
  dev.off()
}

proteionBlast("CALM2")
proteionBlast("CCNI")
proteionBlast("CNOT6L")
proteionBlast("SPIT2")