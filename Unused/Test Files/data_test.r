#TMexampleData TEST

#The purpose of this file is to open TMexampleData and understand how it is structured.

load("/Users/nathanhasegawa/Desktop/NSH_WI26_Rotation/data/TMexampleData.Rdata")

#all.expr: contains the raw gene expression data. Appears to be an array.

#all.meta: has 1115 observations (rows), 10 variables. Columns are:
    #samp: Seems to identify a sample
    #study: likely the study that data point came from (TrTe, V1, V2, V3)
    #cond: ???
    #train: my guess is that it is 1 if the data point was used in a training dataset, 0 otherwise.
    #ID and sID: there are 116 unique entries
    #DLMO25: ???
    #LocalTime: the response variable. It appears to be the time, in hours, when the blood draw
    #was taken (0 is midnight? 23.999 is midnight? These all appear to be rounded, and at first
    #glance I can see overreporting of even numbered miniutes/well-rounded times are overrepresented.)
    #CPdeg: appears to be between 0 and 360. This is probably the angle of the PREDICTED time of day
    #(e.g. 90 corresponds to 6 am, 225 corresponds to 3 pm).
    #CPhrs: it's always CPdeg divided by 15. This makes sense (24 * 15 = 360).
    #THERE ARE 349 NAS IN CPHRS AND CPDEG. I SUSPECT THAT THESE MIGHT BE THE TRAINING DATA.
    #DMLO: Time of Dim Light Melatonin Onset. Essentially a phase marker.

TrTe_data <- all.meta[all.meta$study == "TrTe",]
V1_data <- all.meta[all.meta$study == "V1",]
V2_data <- all.meta[all.meta$study == "V2",]
V3_data <- all.meta[all.meta$study == "V3",]

#cycling.expr: 37 rows and 1115 columns. Each row corresponds to one of the 37 "cycling genes"
#identified in the analysis (e.g. SMAP2, UBE2B, CDC42EP2). These appear to be raw data and not
#studentized.

#Each column corresponds to the gene expression level of each cycling gene in one sample;
#clearly they are represented in log2-space. All the numbers in column 1 appear to be around
#the same few orders of magnitude (7-13); a quick glance at some of the other columns suggest
#that values do not deviate much from this range.

#ratio.expr: contains the data for the ratio of expression of one gene to another. Rows 
#correspond to pairs of genes (e.g. SERTAD1 ZNF101, PHF21A IL13RA1). Columns correspond to
#blood samples. Entries in the columns appear to be roughly between 0 and 3.

#I'm not sure exactly what is contained in ratio.label but it does not appear to contain any
#quantitative entries relevant to the analysis. See how this is called in other files.

#predictor.expr.final: I suspect that this contains the values of the genes in the predictor
#ONLY for each sample. NOT all genes (that's all.expr and all.meta). 
#predictor.expr.nor appears to contain the same data (but I didn't test that in depth; it
#may be false).

#matPredDeg: THIS IS THE HEART OF THE ANALYSIS. The first column is CPdeg (which I think
#is the actual phase of the patient at the time the sample was taken) and the second
#column is matPLpred (which I think is the PREDICTOR-based phase). So the error would be the
#absolute value of their difference (looping around 360).

#From the dimensions of CPtrainZScoreDat and CPtrainRatioDat, it appears that we have 170
#training samples. A quick glance suggests that ALL of the training samples were in the
#TrTe experiment.

#CPtrain contains booleans; these are TRUE if the associated sample is a training sample
#and FALSE otherwise.

#Lots of functions that just convert between different representations of time— use as needed 