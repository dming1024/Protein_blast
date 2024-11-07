BiocManager::install("msa")
library(msa)

system.file("tex", "texshade.sty", package="msa")

#包含有多条蛋白氨基酸序列的fasta文件
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)

#run the clustalW
myFirstAlignment <- msa(mySequences)

#输出pretty的结果
msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
#上述步骤无法成功，需要手动使用tools::texi2pdf命令
tools::texi2pdf("myFirstAlignment.tex",clean = T)
#也会报错，显示fasta not found，需要修改*.tex中的fasta文件，再次试验即可成功

#output为asis在knitr中显示
msaPrettyPrint(myFirstAlignment, y=c(164, 213), output="asis",
               showNames="none", showLogo="none", askForOverwrite=FALSE)


#Phylogenetic Tree
hemoAln2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
library(seqinr)
d <- dist.alignment(hemoAln2, "identity")
library(ape)
hemoTree <- nj(d)
plot(hemoTree, main="Phylogenetic Tree of Hemoglobin Alpha Sequences")


#对alignment结果进行优化
msaPrettyPrint(myFirstAlignment, output="pdf", 
               subset=c(1:6), showNames="left", showLogo="top",
               logoColors="rasmol", shadingMode="similar",
               showLegend=FALSE, askForOverwrite=FALSE)
tools::texi2pdf("myFirstAlignment.tex",clean = T)
#"C:/Users/fan_qiangqiang/AppData/Local/Temp/"

msaPrettyPrint(myFirstAlignment, output="pdf", 
               showNames="left", showLogo="top",
               logoColors="rasmol", shadingMode="similar",
               showLegend=FALSE, askForOverwrite=FALSE)
