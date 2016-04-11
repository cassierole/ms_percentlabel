# This function accepts as argument a csv file (ex."file.csv") containing raw quantitative  
# mass spec data and produces a file in the working directory containing the percentage of  
# C13-labeled phospholipid for each individual species. 
# Experiments carried out on Acinetobacter baumannii fed 2-C13 acetate 
# strains is a character vector of the strains used in the experiment.
percentlabelv3 <- function(datafile,strains,times=c(20,40,60)){
        
        dat <- read.csv(datafile, check.names = FALSE) #import raw mass spec data
        frame <- data.frame() #create an empty data frame for binding
        cname <-names(dat) #store the column names from mass spec data in vector cnames
        rnames <- dat[,grep("name",cname,ignore.case=TRUE)] #vector of sample names
        wt <- grep("WT",rnames,ignore.case=TRUE) #rows containing wt samples

        for (i in 1:(nrow(dat))){
                name <- rnames[i]
                #calculate the percentage labeled for each species and save as variables
                dioPG <- getvalues(i,"795","773",dat,cname)
                paloPG <- getvalues(i,"768","747",dat,cname)
                dipPG <- getvalues(i,"739","719",dat,cname)
                dioPE <- getvalues(i,"764","742",dat,cname)
                paloPE <- getvalues(i,"737","716",dat,cname)
                dipPE <- getvalues(i,"708","688",dat,cname)
                card <- getvalues(i,"711","701",dat,cname)
                
                #create a new data frame containing percentages, then bind this to existing data frame
                newframe <- data.frame("Sample_Name"=name,"PG_C18:1/18:1"=dioPG,"PG_C18:1/16:0"=paloPG,"PG_C16:0/16:0"=dipPG,"PE_C18:1/18:1"=dioPE,"PE_C18:1/16:0"=paloPE,"PE_C16:0/16:0"=dipPE,"Cardiolipin"=card, check.names=FALSE)
                frame <- rbind(frame,newframe) 
        }
        
        membrane <-checkMembrane(rnames)
        frame[,"Membrane"] <- membrane
        time <-checkTime(rnames,times)
        frame[,"Time"] <- time
        ID <-checkID(rnames,strains)
        frame[,"ID"] <-ID
        #creates a file in working directory containing percentages:
        filename <- paste("percent",datafile,sep="_")
        write.csv(frame,file=filename)
        message("Calculation Complete")
        frame
        
}

#a function to calculate percentages
calcper <- function(lab, unlab){
        perlab <- 100 * lab / (lab + unlab)
}

#getvalues accepts as arguments row of database, m/z of parent ions for labeled and unlabeled, and database
#returns percentage of labeled species out of total for that species.
getvalues <-function(row,label,unlabel,data,cnames){
        #store as a vectors the column positions of labeled and unlabeled species
        loclab <- grep(label, cnames, ignore.case=TRUE)
        locunl <- grep(unlabel, cnames, ignore.case=TRUE)
        vallab <- NA
        valunl <- NA
        
        #If there is more than one transition corresponding to a given parent ion, sum the quantitation values
        if ((length(loclab)>0)&&(length(locunl)>0)){
                
                for (i in 1:length(loclab)){
                        vallab <- sum(vallab, data[row,loclab[i]],na.rm=TRUE)  
                }
                
                for (i in 1:length(locunl)){
                        valunl <- sum(valunl, data[row,locunl[i]],na.rm=TRUE)  
                }
                
        }
        calcper(vallab,valunl)
}

#Returns a vector of whether the sample is from the inner or the outer membrane
checkMembrane <- function(rownames){
        memb <- vector(mode = "character", length=length(rownames))
        inner <- grep("_IM",rownames,ignore.case=TRUE)##edited underscores here
        outer <- grep("_OM",rownames,ignore.case=TRUE)##underscores
        for (i in inner){
                memb[i] <- "Inner Membrane"
        }
        for (i in outer){
                memb[i] <- "Outer Membrane"
        }
        memb #If membrane not specified in sample name, entry will be blank. May want to change to NA...
}

#Creates a vector containing the time point of each sample
checkTime <- function(rownames, time){
        ftim <-paste("T",time,sep="") #ex:("T20","T40","T60")
        times <- vector(mode = "character",length = length(rownames))
        for (i in 1:length(ftim)){ #create list, where each object contains a vector of the rows that have the corresponding time.
                timlist <- grep(ftim[i],rownames,ignore.case=TRUE)
                for (j in timlist)
                        times[j] <- paste(time[i],"mins")
        }
        times
}

#Creates a vector which assigns the name of the strain to the sample
#Input the list of sample names and the strains used.
checkID <- function(rownames, strains){
        ID <-vector(mode="character",length=length(rownames))
        for (i in 1:length(strains)){
                posits <- grep(strains[i],rownames,ignore.case=TRUE)
                for (j in posits){
                        ID[j] <- strains[i]
                }
                        
        }
        ID   

}


