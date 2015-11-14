#Last updated: 11/27/07
#Function (in R language) that performs item selection for the development of
#short-forms of scales using an Ant Colony Optmization (ACO) Algorithm.

#Author: Walter Leite, Assitant Professor, Department of Educational Psychology
#University of Florida
#Contact Information: walter.leite@coe.ufl.edu

#The paper "Item selection for the development of short-forms of scales using
#an Ant Colony Optimization algorithm" can be obtained by e-mailing the author.

#This file was created using the Tinn-R 1.19 text editor, available at:
#http://sourceforge.net/projects/tinn-r
#Information about the R language and software download are available at:
#http://www.r-project.org



#Arguments:
#ants: number of ants to send in search of a solution per iteration.
#evaporation: sets the pheromone evaporation rate.
#mplus: name of the MPLUS input file that fits the confirmatory factor analysis model.
#list.items: list containing a vector of items for each factor.
#full:total number of items in the test.
#i.per.f: vector with number of items per factor, in the same order of list.items.
#factors: character vector with names of factors (should be the same in MPLUS input file).
#par: number of parameters estimates with the short-form model
#iter: stopping rule: number of iterations without change.
#resultfile: file where the results for each short form is saved by MPLUS.
#summaryfile: file where the results of each run should be saved by R.
#loc.gammas: vector with the line numbers where the regression coefficients
#of the MIMIC model start and end (locations). Obtain these numbers from the
#TECH1 output of MPLUS.
#loc.variances = vector with the line numbers of the residual variances of the latent factors.
#Also obtain these numbers from the TECH1 output.
#predictors = vector with names of predictors.
#var.predictor = vector with variances of the predictors.

#FUNCTION START:
antcolony = function(ants = 10, evaporation = 0.95, mplus,list.items,full,i.per.f,
  factors,par,iter = 100,resultfile,summaryfile,loc.gammas,loc.variances,
  predictors,var.predictors) {

#creates the table of initial pheromone levels.
include = rep(2,full)

#puts initial best solution (all items selected).
best.so.far.solution = include

#creates a vector with all items.
item.vector = unlist(list.items, use.names = F)

#reads mplus syntax.
setwd("C:/Program Files/Mplus")
input = scan(file = paste(mplus,".inp",sep=""),what="character", sep ="\n", quote = NULL)

#creates a list to store factors.
selected.items = list.items

#starts counting the iterations
count = 1

#starts counting continuous runs regardless of result.
run = 1

#defines initial best so far (overall) pheromone
best.so.far.pheromone = 0
#defines initial best pheromone for the current trial of n ants.
best.pheromone = 0
#defines initial solutions.
previous.solution = include
step = 0

#starts loop through iterations.
while (count <= iter) { 

#sends a number of ants per time.
ant  = 1
while (ant <= ants) {

#selects items for all factors.
for (factor in 1:length(list.items)) {

#selects the items for a short form for the factor
positions = is.element(item.vector,list.items[[factor]])
prob = include[positions]/sum(include[positions])

items = sample(list.items[[factor]], size = i.per.f[factor],replace = F,prob)

#stores selected items.                  
selected.items[[factor]] = items

#replaces the mplus syntax for factor specification.
#make sure the factor name provided is the same as the one in the MPLUS file.
factor.position = grep(paste(factors[factor],"BY"),input,ignore.case=T)
input[factor.position]  = paste(c(factors[factor],"BY",items,";"),collapse =" ")

#finishes loop
}

#creates a vector of selected items.
selected.items = lapply(selected.items,sort)
selected.vector = unlist(selected.items, use.names = F)

#creates a 0/1 vector of the same length of the full form indicating
#whether an item was selected or not for the short form.
select.indicator = is.element(item.vector,selected.vector)
notselect.indicator = (select.indicator == FALSE)

#MODIFY MPLUS SYNTAX.
#when constructing mplus file, leave two lines blank after the usevariables command.

#define categorical variables.
position0 = grep("categorical are",input, ignore.case = T)
input[position0] = paste("categorical are")
input[position0+1] = paste(c(selected.vector[1:round(sum(i.per.f)/2,0)]),
                            collapse = " ")
input[position0+2] = paste(c(selected.vector[(round(sum(i.per.f)/2,0)+1):sum(i.per.f)],";"),
                            collapse = " ")

#make sure there are at least three lines of space between them.                            
#defines variables to be used.                           
position1 = grep("usevariables",input,ignore.case=T)
input[position1] = paste("usevariables are ",predictors)
input[position1+1] = paste(c(selected.vector[1:round(sum(i.per.f)/2,0)]),
                            collapse = " ")
input[position1+2] = paste(c(selected.vector[(round(sum(i.per.f)/2,0)+1):sum(i.per.f)],";"),
                            collapse = " ")

#writes the mplus syntax for fitting the short form.
write.table(input, file = paste(mplus,".inp",sep=""), quote = F,
              col.names = F, row.names = F)

#Run MPLUS with the syntax file specified.
system(paste("mplus.exe ",mplus,".inp",sep=""), wait = TRUE)
outfile = c(paste(mplus,".out",sep=""))

#read the output file to verify the outcome of the analysis.
output = scan(file = outfile, what = character())
#check if the outcome was a non-positive definite matrix.
outcome1 = grep("WARNING",output,ignore.case=T) 
#check if there was non-convergence.
outcome0 = grep("EXCEEDED",output,ignore.case=T)
#set pheromone to zero if the there was non-convergence or a non-positive definite solution.
if ((length(outcome1) == 1) || (length(outcome0) == 1)) {
pheromone = 0 

#writes feedback about non-convergence and non-positive definite.
fit.info = matrix(c(select.indicator,run,count,ant,999,999,round((include),5)),1,)
write.table(fit.info, file = summaryfile, append = T,
              quote = F, sep = " ", row.names = F, col.names = F)

#write feedback about search to a HTML file called iteration.html.
feedback = c(paste("<h1>",run,"-",count,"-",ant,"-",step,"- Failure", "</h1>" ) )
write(feedback, file = "iteration.html")
              
#finishes if for non-convergent cases.
} else { 

#The SAVE command at the MPLUS syntax specifies that results should be
#saved to results.txt.
#scans all results.
results = scan(file=resultfile)
#read the CFI,TLI and RMSEA.
CFI = round(results[(2*par + 4)],3)
TLI = round(results[(2*par + 5)],3)
RMSEA = round(results[(2*par + 7)],3)

#reads the regression coefficients (gammas)
gammas = results[loc.gammas]

#reads the residual variances of the latent variables.
res.variances = results[loc.variances]

#obtains the variance explained by the predictors
var.explained = t(t(gammas)^2*var.predictors)
var.explained = apply(var.explained,1,sum)

#obtains the total variance of the latent variable
#residual variance + variance explained
variances = res.variances + var.explained

#standardizes the gammas.
std.gammas = gammas/sqrt(variances)

#saves information about the selected items and the RMSEA they generated.
fit.info = matrix(c(select.indicator,run,count,ant,CFI,TLI,RMSEA,mean(std.gammas),
  round(include,2)),1,)
 
write.table(fit.info, file = summaryfile, append = T,
              quote = F, sep = " ", row.names = F, col.names = F)
              
#provide feedback about search.
feedback = c(paste("<h1>","run:",run,"count:",count,"ant:",ant,"step:",step,"<br>",
"CFI:",CFI,"TLI:",TLI,"RMSEA:",RMSEA,"<br>",
"GAMMA:",mean(std.gammas), "</h1>" ) )
write(feedback, file = "iteration.html")


#implements fit requirement.
if ((CFI < 0.95)|(TLI <0.95)|(RMSEA > 0.06)) {
pheromone = 0 } else {
#calculate pheromone (mean of standardized gammas).
pheromone = mean(std.gammas)
}

#adjusts count based on outcomes and selects best solution.
if (pheromone >= best.pheromone) {

#updates solution.
best.solution = select.indicator
#updates best RMSEA.
best.RMSEA = RMSEA
#updates best pheromone
best.pheromone = pheromone
} 

#Move to next ant.
ant = ant + 1

#end else clause for converged solutions
}

#ends loop through ants.
}


#adjusts pheromone only if the current pheromone is best than the previous.
if (best.pheromone > best.so.far.pheromone) {

#implements pheromone evaporation.
include = include*evaporation

#Adjusts the pheromone levels.
include.pheromone = best.solution * best.pheromone*(0.1*run)


#updates pheromone.
include = include + include.pheromone

#updates best so far solution and pheromone, with corresponding RMSEA.
best.so.far.solution = best.solution
best.so.far.pheromone = best.pheromone
best.so.far.RMSEA = best.RMSEA
#re-starts count.
count = 1

#end if clause for pheromone adjustment.
} else {

#advances count.
count = count + 1


#adds more pheromone to the best so far solution.
#include.pheromone = best.so.far.solution * best.so.far.pheromone

#updates pheromone.
#include = include + include.pheromone

#finish else clause.
}

#ends loop.
run = run + 1
}

#Compile a matrix witht the final solution.
final.solution = matrix(c(best.so.far.RMSEA,best.so.far.pheromone,best.so.far.solution),1,,
dimnames=list(NULL,c("RMSEA","mean_gamma",item.vector)))

#FUNCTION END.
return(final.solution)
}
#========================================================================

#DESCRIPTION OF FUNCTION OUTPUT:
#The function outputs a vector containing (in this order):
#RMSEA from the model fit of the optimum solution.
#Pheromone level of the optimum solution.
#vector of zeros and ones of the same length of the full form of the scale
#where a 1 indicates that the corresponding item is part of the optimum short form.

#==================================================
#Example:
#Find the optimum 22 item short form from a 39 item scale.

#Arguments:
#ants = 100
#evaporation = 0.90
#mplus = c("d39_short22")
#list.items = list(EM = c("d3","d7","d9","d10","d11","d12","d13","d16","d25",
#        "d29","d32","d33","d34","d35","d36"),
#        DC = c("d1","d4","d5","d14","d15","d17","d18","d24","d27","d28","d31","d39"),
#        AW = c("d2","d6","d8","d22"),
#        SB = c("d19","d20","d26","d37","d38"),
#        SF = c("d21","d23","d30"))
#full = 39
#i.per.f = c(5,5,4,5,3)
#factors = c("EM","DC","AW","SB","SF")
#par = 47
#resultfile = c("results_short_form22.txt")
#summaryfile = c("output.txt")
#loc.gammas = c(18)
#loc.variances = c(33,35,38,42,47)
#predictors = c("pcsugar")
#var.predictors = 0.138


#run the function
#short22 = antcolony(ants,evaporation,mplus,list.items,full,i.per.f,factors,par,
#                      steps,resultfile,summaryfile,loc.gammas,loc.variances,
#                      predictors,var.predictors)



