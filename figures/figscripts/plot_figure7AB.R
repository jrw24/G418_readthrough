### Rscript for using Xtail 

library("xtail")
library("glue")

args <- commandArgs(trailingOnly = TRUE)

libSetPathForR <- args[1]
rootDir <- args[2]



samplelist <- c(
	"1_allAG_1_Untr_A",
	"2_allAG_2_Untr_B",
	"3_allAG_3_G418_500_A",
	"4_allAG_4_G418_500_B",
	"5_allAG_5_G418_2000_A",
	"6_allAG_6_G418_2000_B",
	"7_allAG_17_G418_10min_A",
	"8_allAG_18_G418_10min_B"
	)

samplelist_renamed <- gsub("(.*?)_(.*)", "\\2_\\1", samplelist)
print(samplelist_renamed)




mrnaTables <- glue("{rootDir}/Data/RNA/FPassignment/hg38_protCode/HEK_G418/analysis/countTablesRAW/HEK_G418_mRNA_RAWcounts.csv")
rpfTables <- glue("{rootDir}/Data/RPF/FPassignment/hg38_protCode/allG418/analysis/countTablesRAW/allG418_cds_RAWcounts.csv")
UTRfilestring <- glue("{rootDir}/genomes/gencodeV30_protCode_TermStopCodon_validUTRs_UTRs.csv")


outDir <- glue("{rootDir}/Data/RNA/FPassignment/hg38_protCode/HEK_G418/analysis/Xtail")
dir.create(outDir, showWarnings=FALSE)


mrna <- read.table(mrnaTables, sep=',', stringsAsFactors = FALSE, header = TRUE)
colnames(mrna) <- c('tr_id',samplelist_renamed)
rpf <- read.table(rpfTables, sep=',', stringsAsFactors = FALSE, header = TRUE)

# rpf 

rpf <- rpf[,c(1:7, 18, 19)] ## first column in now 'tr_id'
head(mrna)
head(rpf)
colnames(rpf) <- c('tr_id', samplelist_renamed)


conditions <- c("control", 'g418_500', 'g418_2K', 'g418_10min')
condition <- rep(conditions, each=2)
condition



#### --> G418 2K
cntrl <- "control"
treat <- "g418_2K"

mrnaY <- mrna[,c(1, 2:3, 6,7)]
mrnaX <- mrnaY[,-1]
rownames(mrnaX) <- mrnaY$tr_id
rpfY <- rpf[,c(1, 2:3, 6,7)]
rpfX <- rpfY[,-1]
rownames(rpfX) <- rpfY$tr_id
condX <- condition[c(1:2, 5:6)]

head(mrnaX)
head(rpfX)
condX

# test.results <- xtail(mrnaX, rpfX, condX, bins=1000)
test.results <- xtail(mrna = mrnaX, 
                       rpf = rpfX, 
                       condition = condX, 
                       minMeanCount = 50,
                       ci = 0.95,
                       normalize = TRUE,
                       threads = 40,
                       bins=10000,
                      )



### extract results table:
resTab <- resultsTable(test.results, log2FCs = TRUE, log2Rs = TRUE)
head(resTab, 5)
UTRtab <- read.csv(UTRfilestring, header=TRUE, stringsAsFactors = FALSE)
resTab$tr_id <- rownames(resTab)

ptab <- merge(resTab, UTRtab, by.x = 'tr_id', by.y = 'X.transcript')
head(ptab)

out1 <- glue('{outDir}/Untr_vs_g4182K_XtailResTab.csv')
write.csv(ptab, file=out1, row.names=TRUE)


#### --> G418 10 min
cntrl <- "control"
treat <- "g418_10min"

mrnaY <- mrna[,c(1, 2:3, 8:9)]
mrnaX <- mrnaY[,-1]
rownames(mrnaX) <- mrnaY$tr_id
rpfY <- rpf[,c(1, 2:3, 8:9)]
rpfX <- rpfY[,-1]
rownames(rpfX) <- rpfY$tr_id
condX <- condition[c(1:2, 7:8)]


test.results2 <- xtail(mrna = mrnaX, 
                       rpf = rpfX, 
                       condition = condX, 
                       minMeanCount = 50,
                       ci = 0.95,
                       normalize = TRUE,
                       threads = 40,
                       bins=10000,
                      )

resTab2 <- resultsTable(test.results2, log2FCs = TRUE, log2Rs = TRUE)
UTRtab <- read.csv(UTRfilestring, header=TRUE, stringsAsFactors = FALSE)
resTab2$tr_id <- rownames(resTab2)

# head(resTab2)

ptab2 <- merge(resTab2, UTRtab, by.x = 'tr_id', by.y = 'X.transcript')
out2 <- glue('{outDir}/Untr_vs_g41810min_XtailResTab.csv')
write.csv(ptab2, file=out2, row.names=TRUE)


### --> G418 500
cntrl <- "control"
treat <- "g418_500"

condition

mrnaY <- mrna[,c(1, 2:3, 4,5)]
mrnaX <- mrnaY[,-1]
rownames(mrnaX) <- mrnaY$tr_id
rpfY <- rpf[,c(1, 2:3, 4,5)]
rpfX <- rpfY[,-1]
rownames(rpfX) <- rpfY$tr_id
condX <- condition[c(1:2, 3:4)]


# test.results <- xtail(mrnaX, rpfX, condX, bins=1000)
test.results3 <- xtail(mrna = mrnaX, 
                       rpf = rpfX, 
                       condition = condX, 
                       minMeanCount = 50,
                       ci = 0.95,
                       normalize = TRUE,
                       threads = 40,
                       bins=10000,
                      )

resTab3 <- resultsTable(test.results3, log2FCs = TRUE, log2Rs = TRUE)
UTRtab <- read.csv(UTRfilestring, header=TRUE, stringsAsFactors = FALSE)
resTab3$tr_id <- rownames(resTab3)


ptab3 <- merge(resTab3, UTRtab, by.x = 'tr_id', by.y = 'X.transcript')

out3 <- glue('{outDir}/Untr_vs_g418_500_XtailResTab.csv')
write.csv(ptab3, file=out3, row.names=TRUE)





