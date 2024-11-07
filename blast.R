
#multiple blast of proteins
#包含有多条蛋白氨基酸序列的fasta文件
library(msa)
mySequenceFile <- "./ncbi/CREBBP.fasta"
mySequences <- readAAStringSet(mySequenceFile)

#run the clustalW
myFirstAlignment <- msa(mySequences)
#输出pretty的结果
msaPrettyPrint(myFirstAlignment, output="pdf", file="result.pdf",
               showNames="left", shadingMode="identical",
               showLogo = "none",consensusColor="ColdHot",
               #shadingModeArg = "standard area",
               furtherCode = c("\\showruler{1}{top}"),
               showLegend=FALSE, askForOverwrite=FALSE)

tools::texi2pdf("result.tex",clean = T)
#"C:/Users/fan_qiangqiang/AppData/Local/Temp/"

#输出HAT domain区域
msaPrettyPrint(myFirstAlignment, output="pdf", file="resultHAT.pdf",
               showNames="left", showLogo="top",y=c(1300, 1670),
               logoColors="rasmol", shadingMode="identical",
               #shadingModeArg = "standard area",
               showLegend=FALSE, askForOverwrite=FALSE)
tools::texi2pdf("resultHAT.tex",clean = T)
#
#Phylogenetic Tree
hemoAln2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
library(seqinr)
d <- dist.alignment(hemoAln2, "identity")
library(ape)
hemoTree <- nj(d)
pdf("phylogenetic_of_CREBBP.pdf")
plot(hemoTree, main="Phylogenetic Tree of CREBBP")
dev.off()

#ep300基因分析
########
mySequenceFile_EP300 <- "./ncbi/EP300.fasta"
mySequences_EP300 <- readAAStringSet(mySequenceFile_EP300)

#run the clustalW
myFirstAlignment_EP300 <- msa(mySequences_EP300)
#输出pretty的结果
msaPrettyPrint(myFirstAlignment_EP300, output="pdf", file="result_EP300.pdf",
               showNames="left", shadingMode="identical",
               showLogo = "none",consensusColor="ColdHot",
               #shadingModeArg = "standard area",
               furtherCode = c("\\showruler{1}{top}"),
               showLegend=FALSE, askForOverwrite=FALSE)

tools::texi2pdf("result_EP300.tex",clean = T)
#"C:/Users/fan_qiangqiang/AppData/Local/Temp/"

#输出HAT domain区域
msaPrettyPrint(myFirstAlignment_EP300, output="pdf", file="resultHAT.pdf",
               showNames="left", showLogo="top",y=c(1280, 1600),
               logoColors="rasmol", shadingMode="identical",
               #shadingModeArg = "standard area",
               showLegend=FALSE, askForOverwrite=FALSE)
tools::texi2pdf("resultHAT.tex",clean = T)
#
#Phylogenetic Tree
hemoAln2 <- msaConvert(myFirstAlignment_EP300, type="seqinr::alignment")
library(seqinr)
d <- dist.alignment(hemoAln2, "identity")
library(ape)
hemoTree <- nj(d)
pdf("phylogenetic_of_EP300.pdf")
plot(hemoTree, main="Phylogenetic Tree of EP300")
dev.off()

#将4个物种的2个基因，放在一个里面进行blast
mySequenceFile_both <- "./ncbi/both.fasta"
mySequences_both <- readAAStringSet(mySequenceFile_both)

#run the clustalW
myFirstAlignment_both<- msa(mySequences_both)
#输出pretty的结果
msaPrettyPrint(myFirstAlignment_both, output="pdf", file="result_both.pdf",
               showNames="left", shadingMode="identical",
               showLogo = "none",consensusColor="ColdHot",
               #shadingModeArg = "standard area",
               furtherCode = c("\\showruler{1}{top}"),
               showLegend=FALSE, askForOverwrite=FALSE)

tools::texi2pdf("result_both.tex",clean = T)
#"C:/Users/fan_qiangqiang/AppData/Local/Temp/"

#输出HAT domain区域
msaPrettyPrint(myFirstAlignment_both, output="pdf", file="resultHAT_both.pdf",
               showNames="left", showLogo="top",y=c(1280, 1700),
               logoColors="rasmol", shadingMode="identical",
               #shadingModeArg = "standard area",
               showLegend=FALSE, askForOverwrite=FALSE)
tools::texi2pdf("resultHAT_both.tex",clean = T)
#
#Phylogenetic Tree
hemoAln2 <- msaConvert(myFirstAlignment_both, type="seqinr::alignment")
library(seqinr)
d <- dist.alignment(hemoAln2, "identity")
library(ape)
hemoTree <- nj(d)
pdf("phylogenetic_of_both.pdf")
plot(hemoTree, main="Phylogenetic Tree of Both Genes")
dev.off()


#MSH分析
#将4个物种的2个基因，放在一个里面进行blast
mySequenceFile_both <- "msh.fasta"
mySequences_both <- readAAStringSet(mySequenceFile_both)

#run the clustalW
myFirstAlignment_both<- msa(mySequences_both)
#输出pretty的结果
msaPrettyPrint(myFirstAlignment_both, output="pdf", file="result_msh.pdf",
               showNames="left", shadingMode="identical",
               showLogo = "none",consensusColor="ColdHot",
               #shadingModeArg = "standard area",
               furtherCode = c("\\showruler{1}{top}"),
               showLegend=FALSE, askForOverwrite=FALSE)

tools::texi2pdf("result_msh.tex",clean = T)
#"C:/Users/fan_qiangqiang/AppData/Local/Temp/"

#输出MSH domain区域
msaPrettyPrint(myFirstAlignment_both, output="pdf", file="resultHAT_msh.pdf",
               showNames="left", showLogo="top",y=c(500, 700),
               logoColors="rasmol", shadingMode="identical",
               #shadingModeArg = "standard area",
               showLegend=FALSE, askForOverwrite=FALSE)
tools::texi2pdf("resultHAT_msh.tex",clean = T)
#
#Phylogenetic Tree
hemoAln2 <- msaConvert(myFirstAlignment_both, type="seqinr::alignment")
library(seqinr)
d <- dist.alignment(hemoAln2, "identity")
library(ape)
hemoTree <- nj(d)
pdf("phylogenetic_of_msh.pdf")
plot(hemoTree, main="Phylogenetic Tree of MSH6")
dev.off()
