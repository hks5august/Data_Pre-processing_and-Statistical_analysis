#################################Survival Analysis  on the basis of expression of TCGA #############
library(survival)
install.packages("rpart")
# install.packages('survminer')
install.packages("ranger")
library(ranger)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
#read matrix#####
matrix_os <- read.csv('/Users/anjali/Harpreet_project/survival/cancer_dfi.csv', header = T) ### input_matrix####
#data("ovarian")
dim(matrix_os)
data<-as.data.frame(matrix_os[,2:29])
###Drop NA ##################
data1 <- na.omit(data)
dim(data1)

glimpse(data1) #####all_columns
#glimpse(ovarian)
#help("ovarian")
#data$D_PFS_FLAG <- factor(data$D_PFS_FLAG, levels = c("alive", "dead"), labels = c("0", "1"))

######## fit survival data using Kaplen-Meier method #############
#####################Overall survival ##############################
surv_os <- Surv(time = data1$OS.time, event = data1$OS) 
##################Disease specific survival ##################
surv_DSS <- Surv(time = data1$DSS.time, event = data1$DSS) 
##################Progression free survival #########################
surv_PFI <- Surv(time = data1$PFI.time, event = data1$PFI) 
##################Disease free survival #########################
surv_DFI<- Surv(time = data1$DFI.time, event = data1$DFI)


############### Kaplan plot ##########################



################### Mean  #########################

for (i in colnames(data1[3:28]))
{
  
  jj <- survfit(surv_PFI ~ (data1[[i]])>mean((data1[[i]])),data = data1) #########Fit survival curve for kaplan plot ########
  cox <- coxph(surv_PFI ~ (data1[[i]])>mean((data1[[i]])), data = data1)   ###########calculation of cox value and HR###################
#############   OUTPUT IN THE FORM OF HR AND P-VALUE ################  
  write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
              file="survival/PFI_greatermean.csv",row.names=F,col.names=F,sep = ',',append = T)
############ OUTPUT FOR SIGNIFICANT GENES #####################################################################  
  if ((coef(summary(cox))[5]<0.05) & (!is.na(coef(summary(cox))[5])) & (!is.na(coef(summary(cox))[2])))
  {write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
               file="survival/PFI_greatermean_significant.csv",row.names=F,col.names=F,sep = ',',append = T)}
  
}
for (i in colnames(data1[3:28]))
{
  
  jj <- survfit(surv_PFI ~ (data1[[i]])<mean((data1[[i]])),data = data1)
  cox <- coxph(surv_PFI ~ (data1[[i]])<mean((data1[[i]])), data = data1)
  
  write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
              file="survival/PFI_lessmean.csv",row.names=F,col.names=F,sep = ',',append = T)
  if ((coef(summary(cox))[5]<0.05) & (!is.na(coef(summary(cox))[5])) & (!is.na(coef(summary(cox))[2])))
  {write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
               file="survival/PFI_lessmean_significant.csv",row.names=F,col.names=F,sep = ',',append = T)}
  
}
############## Median ##################################
for (i in colnames(data1[3:28]))
{
  
  jj <- survfit(surv_PFI ~ (data1[[i]])<median((data1[[i]])),data = data1) #########Fit survival curve for kaplan plot ########
  cox <- coxph(surv_PFI  ~ (data1[[i]])<median((data1[[i]])), data = data1) ###########calculation of cox value and HR###################
  #############   OUTPUT IN THE FORM OF HR AND P-VALUE ################  
  write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
              file="survival/PFI_comp_lessmedian.csv",row.names=F,col.names=F,sep = ',',append = T)
  ############ OUTPUT FOR SIGNIFICANT GENES #####################################################################  
  if ((coef(summary(cox))[5]<0.05) & (!is.na(coef(summary(cox))[5])) & (!is.na(coef(summary(cox))[2])))
  {write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
               file="survival/PFI_comp_lessmedian_significant.csv",row.names=F,col.names=F,sep = ',',append = T)}
  
}
for (i in colnames(data1[3:28]))
{
  
  jj <- survfit(surv_PFI ~ (data1[[i]])>median((data1[[i]])),data = data1)
  cox <- coxph(surv_PFI ~ (data1[[i]])>median((data1[[i]])), data = data1)
  
  write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
              file="survival/PFI_comp_GREATERmedian.csv",row.names=F,col.names=F,sep = ',',append = T)
  if ((coef(summary(cox))[5]<0.05) & (!is.na(coef(summary(cox))[5])) & (!is.na(coef(summary(cox))[2])))
  {write.table(cbind(i,jj[[1]][1],jj[[1]][2],coef(summary(cox))[2],coef(summary(cox))[5]), 
               file="survival/PFI_comp_GREATERmedian_significant.csv",row.names=F,col.names=F,sep = ',',append = T)}
  
}


####################for single gene######################

jj <- survfit(surv_DSS ~ data1$ > mean(data1$DTL),data = data1) ###Fit model ######

cox <- coxph(surv_DSS ~ data1$DTL  > mean(data1$DTL), data = data1) ########COX value in the form of HR ################
################# SUMMARY OF COX values #################################################

summary(cox)
jj

###############SURVIVAL PLOTS##############################################################################
##################ggsurvplot################################################################################
pp_plot<- ggsurvplot(jj, legend.title = "Expression",data = data1,title="                                     DSS curve for DTL",xlab = "Time in days", legend=c(0.8,0.8), risk.table = TRUE,break.time.by = 500,legend.labs = c(" >mean", " < mean"),risk.table.y.text.col = T, risk.table.y.text = FALSE)
pp$plot <- pp$plot+ 
  ggplot2::annotate("text", x = 1000, y = 0.15, # x and y coordinates of the text
                    label = "HR =2.255\n p-value = 0.000467", size = 4)

pp_plot
