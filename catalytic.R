# load dataset

master_data <- read.csv("C:\\Users\\atsui\\Dropbox\\Adrian\\Cambridge\\Part II Genetics\\Project\\GeocodedDatabaseDengueQSNICH.csv")

# install packages

install.packages("lubridate")
library("lubridate")

# cleanup -----------------------------------------------------------------

# cleanup into "age of adm" | "year of adm" | "serostatus"

  # numerise Interpretation: 
    # susceptibles classed as 0; 
    # primary infections as 1; 
    # secondary infections as 2; 
    # all others discarded

    master_data_Classified <- master_data
    master_data_Classified$Interpretation[master_data_Classified$Interpretation == "No Evidence of Recent Flavivirus Infection"] <- 0
    master_data_Classified$Interpretation[master_data_Classified$Interpretation == "Acute Primary Dengue Infection"] <- 1
    master_data_Classified$Interpretation[master_data_Classified$Interpretation == "Acute Secondary Dengue Infection"] <- 2
    master_data_Classified <- subset(master_data_Classified, 
                                     master_data_Classified$Interpretation == 0 | 
                                      master_data_Classified$Interpretation == 1 | 
                                       master_data_Classified$Interpretation == 2)

  # compute ages at admission

    # convert date strings into dates
    master_data_Classified$Birth.date <- dmy(master_data_Classified$Birth.date)
    master_data_Classified$Admission.Date <- dmy(master_data_Classified$Admission.Date)
    
    # compute and add age-of-admission column
    age_adm <- floor(time_length(difftime(master_data_Classified$Admission.Date, master_data_Classified$Birth.date), "years"))
    master_data_Classified$age_adm <- age_adm
    
  # final cleanup
    master_data_Classified <- data.frame(master_data_Classified$age_adm, 
                                         master_data_Classified$Admission.Date,
                                         master_data_Classified$Interpretation)
    colnames(master_data_Classified) <- c("age_adm", "year_adm", "serostatus")
    master_data_Classified$year_adm <- as.numeric(format(master_data_Classified$year_adm, "%Y"))
    master_data_Classified$serostatus <- as.numeric(master_data_Classified$serostatus)

# seroprevalence, constant foi --------------------------------------------

# compute seroprevalence
    
    # 1. assuming constant FOI through time
        
      # fit GLM
        
        # conglomerate all seropositives
        
        master_data_Classified_Simple <- master_data_Classified
        master_data_Classified_Simple$serostatus[master_data_Classified_Simple$serostatus == 2] <- 1
        master_data_Classified_Simple <- master_data_Classified_Simple[master_data_Classified_Simple$age_adm != 0,]
      
        CAT_constantFOI <- glm(formula = serostatus ~ 1 + log(age_adm), 
                               data = master_data_Classified_Simple,
                               family = binomial(link = "cloglog")
                               )
        
        summary(CAT_constantFOI)
        
      # overlay data onto catalytic curve
      
        # initialise output matrix
        
        SP = matrix(0,
                    ncol = 2,
                    nrow = max(master_data_Classified$age_adm) - min(master_data_Classified$age_adm) + 1)
        
        SP[,2] <- min(master_data_Classified$age_adm):max(master_data_Classified$age_adm)
        
        # compute seroprevalence for ages
        
        for (a in min(master_data_Classified$age_adm):max(master_data_Classified$age_adm)){
          
          SP[a,1] <- 1 - nrow(master_data_Classified[master_data_Classified$age_adm == a 
                                                     & master_data_Classified$serostatus == 0,])/
            nrow(master_data_Classified[master_data_Classified$age_adm == a,])
          
        }
        
        SP <- as.data.frame(SP)
        colnames(SP) <- c("seroprevalence", "age")
        
        # final plot
        
        ages = 1:20
        plot(ages, 1 - exp(-0.44774*ages), type = "l", xlab = "age / years", ylab = "seroprevalence")
        points(SP$age, SP$seroprevalence)