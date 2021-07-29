# Tumor_Microenviornment

# 1.Differential expression analysis

## 1.1 Scirpt
DEA.R
## 1.2 input
① IDs of control gorup (.txt file) <br>
② IDs of treat group (.txt file) <br>
③ TCGA project name (the script would download both clinical and expression matrix automatically with the project name)<br>
<br>example:<br>
` 
Rscript DEA.R ./test/control.txt ./test/treat.txt TCGA-ESCA
` 
## 1.3 output
output files listed below: <br>
>result
>>DEA <br>
>>>DEG.txt: differential expression analysis result file <br>
>>>heatmap.png: heatmap for DEGs <br>
>>>volcano.png: volcano for DEGs <br>

>result<br>
>>Stemness<br>
>>>stemness.txt: stemness score of each patients <br>

# 2.TME
## 2.1 Scirpt
Immune.R
## 2.2 Input
① RNA-Seq file path
② IDs of control gorup (.txt file) <br>
③ IDs of treat group (.txt file) <br>
④ parameter：mean OR median to define score cutoff <br>
⑤ survival parameter: OS/DFI/DSS/PFI <br>
⑥ gene names that you focus <br>
 <br>example: <br>
`
Rscript.exe Immune.R ./result/RNAExp.txt ./test/control.txt ./test/treat.txt mean OS SOX2,TP63
`
## 2.3 Output
>result
>>Immune
>>>estimate_score.txt:Stromal and Immune Score of each patient <br>
>>>estimate-ttest.png:violin plot of Stromal and Immune score of two groups (input groups) with p-value(student's t test) 
>>>estimate-survival.png: two groups (divided by mean/median) survival analysis with p-value<br>
>>>xCell_ttest.txt: student's t test result of each cell compounds between input groups. Four columns including cell name, mean of control group, mean of treat group and p-value<br>

>result<br>
>>xCell-ttest: visualizaiton of t-test result which p-value less than 0.05 <br>
>>Gene name (example:SOX2): 64 pngs of correlation plots between gene expression and cell compounds.
