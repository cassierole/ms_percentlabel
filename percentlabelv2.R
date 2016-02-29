# This function accepts as argument a csv file (ex."file.csv") containing raw quantitative  
# mass spec data and produces a file in the working directory containing the percentage of  
# C13-labeled phospholipid for each individual species. 
# Experiments carried out on Acinetobacter baumannii fed 2-C13 acetate 
percentlabelv2 <- function(datafile){
        
        dat <- read.csv(datafile, check.names = FALSE) #import raw mass spec data
        frame <- data.frame() #create an empty data frame for binding
        cname <-names(dat) #store the column names from mass spec data in vector cnames
        
        
        for (i in 1:(nrow(dat))){
                name <- dat[i,grep("name",cname,ignore.case=TRUE)]
                #calculate the percentage labeled for each species and save as variables
                dioPG <- getvalues(i,"795","773",dat,cname)
                paloPG <- getvalues(i,"768","747",dat,cname)
                dipPG <- getvalues(i,"739","719",dat,cname)
                dioPE <- getvalues(i,"764","742",dat,cname)
                paloPE <- getvalues(i,"737","716",dat,cname)
                dipPE <- getvalues(i,"708","688",dat,cname)
                card <- getvalues(i,"711","701",dat,cname)

                #create a new data frame containing percentages, then bind this to existing data frame
                newframe <- data.frame("Sample_Name"=name,"PG C18:1/18:1"=dioPG,"PG C18:1/16:0"=paloPG,"PG C16:0/16:0"=dipPG,"PE C18:1/18:1"=dioPE,"PE C18:1/16:0"=paloPE,"PE C16:0/16:0"=dipPE,"Cardiolipin"=card, check.names=FALSE)
                frame <- rbind(frame,newframe) 
        }
        
        #creates a file in working directory containing percentages:
        filename <- paste("percent_labeled",datafile,sep="_")
        write.csv(frame,file=filename)
        print("Calculation Complete")
        #frame
        
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
        