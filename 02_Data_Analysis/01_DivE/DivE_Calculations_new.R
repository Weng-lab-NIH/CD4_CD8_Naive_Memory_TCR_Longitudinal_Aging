setwd('C:/Users/super/Desktop/DivE_Modified/')
suppressWarnings(suppressMessages(library(DivE)))
suppressWarnings(suppressMessages(library(argparse)))

#DESCRIPTION OF SCRIPT

docstring <- "DESCRIPTION: This R script utilizes the DivE package (v1.0) in order to estimate the \\n\\
number of unique TCRs for a given person. \\n\\
This script is broken down into the following steps: \\n\\n\\
1) Generate the rarefaction curve based on the normalized rarefaction length (lowest UMI count for 1st visit alpha/beta & 2nd visit alpha/beta) \\n\\
2) Generate the rarefaction curve based on half of the normalized rarefaction length \\n\\
3) Fit the specified model to the rarefaction curve \\n\\
4) Project the fit model to the specified cell count (100x less than actual cell count) \\n\\
5) Use the fit model data to calculate to the actual cell count \\n\\
6) Use the fit model data to calculate the normalized cell count \\n\\
7) Save the variables generated for future use into an .Rdata file \\n\\
"

#Create Parser Object
parser <- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')
parser$add_argument("input_directory",help = "Specify the current working directory",type = "character")
parser$add_argument("input_file",help = "Specify the DivE formatted file that you want to use for calculations",type = "character")
parser$add_argument("sample_name",help = "Specify the sample name associated with this file",type = "character")
parser$add_argument("rarefaction_iteration",
                    help = "Specify the UMI number that DivE will iterate to the length of the sample to generate the rarefaction length",
                    type = 'integer')
parser$add_argument("max_rarefaction_length",
                    help = "Specify the maximum UMI length that DivE will calculate to for this sample",
                    type = 'integer')
parser$add_argument("cell_count",help = "Specify the actual cell count of the sample",type = "integer")

parser$add_argument("normalized_cell_count",help = "Specify the normalized cell count for subsets or CD4-CD8",
                    type = 'integer')
parser$add_argument("model_number",
                    help = "Specify the model number that you want to use for calculations",type = 'integer')
args <- parser$parse_args()

#Initialize model data & intiate time calculations
start_timeperfile <-Sys.time()
sample <- read.delim(args$input_file,sep = '\t',stringsAsFactors = FALSE,header = FALSE)
data(ModelSet)
data(ParamRanges)
data(ParamSeeds)
testmodels <- list()
testmeta <- list()
paramranges <- list()

#Choose specific model data
model_name<-names(ModelSet[args$model_number])
testmodels <- ModelSet[args$model_number]
testmeta <- ParamSeeds[args$model_number]
paramranges<- ParamRanges[args$model_number]


#Create subsamples
dss_1 <- divsubsamples(sample, nrf=args$rarefaction_iteration, minrarefac=1, 
                       maxrarefac=args$max_rarefaction_length, NResamples=1000)
dss_2 <- divsubsamples(sample, nrf=args$rarefaction_iteration, minrarefac=1, 
                       maxrarefac=args$max_rarefaction_length*.5, NResamples=1000)
dss <- list(dss_1,dss_2)


#Generate Rarefaction Data, Fit Model, & Project to 1/100 Cell Count
first_cell_proj <- round(args$cell_count/100)
result <- DiveMaster(models=testmodels, init.params=testmeta, param.ranges = paramranges,
                     main.samp=sample, subsizes=c(args$max_rarefaction_length, args$max_rarefaction_length*.5),dssamps = dss,
                     NResamples=1000, nrf=args$rarefaction_iteration, fitloop=2, numit=1e3,tot.pop = first_cell_proj)
subsample_estimated_div <- result$estimate[1]
model_score <- result$model.scores[1]

#Check if model fit is plausible
pc1 <- result$fmm[[model_name]]$plausibility[1,2]
pc2 <- result$fmm[[model_name]]$plausibility[2,2]

if (pc1 + pc2 == 2){
  model_fit_plausible <- TRUE
}else{
  model_fit_plausible <- FALSE
  }

#Calculate curvature of subsamples
curv1 <- Curvature(dss_1)
curv2 <- Curvature(dss_2)


#Project to actual cell count & normalized cell count
actual_estimated_div <- popdiversity(result,popsize = args$cell_count)$estimate
norm_estimated_div <- popdiversity(result,popsize = args$normalized_cell_count)$estimate


#Check if actual projection compared to 100x projection makes sense
if (actual_estimated_div > subsample_estimated_div){
  actual_sub_projection <- TRUE
}else{
  actual_sub_projection <- FALSE
}

#Check if actual projection makes sense compared to normalized cell projection
if (args$cell_count >= args$normalized_cell_count & actual_estimated_div >= norm_estimated_div){
  actual_norm_projection <- TRUE
}else if (args$cell_count <= args$normalized_cell_count & actual_estimated_div <= norm_estimated_div){
  actual_norm_projection <- TRUE
}else{
  actual_norm_projection <- FALSE
}

#Check if 100x projection makes sense compared to normalized cell projection
if ((args$cell_count/100) >= args$normalized_cell_count & subsample_estimated_div >= norm_estimated_div){
  sub_norm_projection <- TRUE
}else if ((args$cell_count/100) <= args$normalized_cell_count & subsample_estimated_div <= norm_estimated_div){
  sub_norm_projection <- TRUE
}else{
  sub_norm_projection <- FALSE
}



print(result) # Comparison of combined scores
print(summary(result)) # Summary of combined score
print(result$ssm) # Detailed comparison of scores
print(result$fmm) #summary of fitting process to model

#Draw Rarefaction Curves
image_1 <- paste(args$sample_name,args$model_number,"sample_range-plot.png",sep='-')
image_1 <- file.path(getwd(),"Plots",image_1)
png(file = image_1)
plot(result$fmm[[model_name]])
dev.off()

image_2 <- paste(args$sample_name,args$model_number,"global_range-plot.png",sep='-')
image_2 <- file.path(getwd(),"Plots",image_2)
png(file = image_2)
plot(result$fmm[[model_name]], range = 'global')
dev.off()


end_timeperfile <- Sys.time()
time_takenperfile <- difftime(end_timeperfile,start_timeperfile,units = 'mins')
print(time_takenperfile)

#Save Results
results_df <- data.frame(sample.name = args$sample_name,model.number = args$model_number,model.name = model_name,
                         bootstrap.length = 100,rarefaction.iteration = args$rarefaction_iteration,
                         max.rarefaction.length = args$max_rarefaction_length,subsampled.100x.cell.pcount = first_cell_proj,
                         actual.cell.count = args$cell_count,normalized.cell.count = args$normalized_cell_count,
                         model.score = model_score,subsampled.100x.estimation = subsample_estimated_div,
                         actual.estimation = actual_estimated_div,normalized.estimation = norm_estimated_div,curvature.main.sample = curv1,
                         curvature.half.sample = curv2,average.curvature = mean(c(curv1,curv2)),
                         model.fit.check = model_fit_plausible,actual.subsample.check = actual_sub_projection,
                         actual.normalized.projection = actual_norm_projection,subsample.normalized.projection = sub_norm_projection,
                         minutes.to.calculations = time_takenperfile)
summary_name <- paste(args$sample_name,"DivE_summary.csv",sep = '-')
summary_name <- file.path(getwd(),"DivE_Summaries",summary_name)
write.csv(results_df,summary_name)
r_object <- paste(args$sample_name,"DivE_environment.Rdata",sep = '-')
r_object <- file.path(getwd(),"R_objects",r_object)
save.image(r_object)
