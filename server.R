
source("LoadHelperVariance.R") 
    
library(shiny)        

  library(knitr)     
library(topGO)  
library(dplyr)     
library(DT)  
library(Hmisc)  
library("shinyURL")     
 
 phenolevels=c("Insufficient data","L. essential","Sig slow","L. dispensable","Sig fast","Unselected")  
 phenolevelslong=c("Insufficient data","Likely essential","Significantly slow","Likely dispensable","Significantly fast","Unselected")  
  phenolevelscolor=c("black", "#ed1c24", "blue","#007603","purple", "darkgray")
 names(phenolevelscolor)=phenolevels
     names(phenolevelslong)=phenolevels    
	 fast=F
	   if(fast==FALSE){
files <- list.files(path = "./countfiles/", pattern = "*.csv", full.names = T, recursive = FALSE)


inputs <- lapply(files, loadSTM)  #Iterate through files loading each to make a list of dataframes
tables <- vector("list", length(inputs))
arrays <- vector("list", length(inputs))
for (i in 1:length(inputs)) {
    tables[[i]] <- inputs[[i]]$table
    arrays[[i]] <- inputs[[i]]$extra
}
fullSet <- ldply(tables, data.frame)  #Merge these dataframes into one massive one
save(file="cache",list=c("tables","arrays","fullSet","files"))
}else{
load(file="cache")
}
fullSet$Relative.Growth.Rate=fullSet$fitness

 
experimentTable<-as.data.frame(table(fullSet$experiment))
experimentTable<-data.frame(Experiment=experimentTable$Var1,Genes=experimentTable$Freq)
 


otherdata = data.frame(thingie = -5:10)
otherdata$halfwidth = 2 * sqrt(exp(-otherdata$thingie))
otherdata$xmin = 1.2 - otherdata$halfwidth
otherdata$xmax = 1.2 + otherdata$halfwidth
# ggplot(comb[comb$Consensus.phenotype!='',],aes(y=-log(variance),x=fitness,color=barseq.phenotype))+geom_point(alpha=1)+coord_cartesian(xlim=c(0,1.5),ylim=c(-5,10))+geom_errorbarh(data=otherdata,aes(xmin=xmin,xmax=xmax,y=thingie,x=xmin+xmax/2),color='gray')+scale_color_manual(values=c('white','blue','red','green','yellow'))
geneinfo <- read.csv("./otherdata/geneinfo.csv")
genomiclocations <- read.csv("./otherdata/genomiclocations.csv", header=TRUE)
genomiclocations$csome=as.numeric(as.character(genomiclocations$csome))

cloneIDs <- read.csv("./otherdata/stm_analyses_repository_copy_150616_2.csv", header=TRUE)
barseqtemp <- read.table("./otherdata/barseqlookup.txt",sep="\t", header=TRUE)

homology<- read.csv("./otherdata/VectorDataGC", header=TRUE)
 

addExtraData<-function(df){
fullSet2<-df
fullSet2$lower=fullSet2$Relative.Growth.Rate-2*sqrt(fullSet2$variance)
	fullSet2$upper=fullSet2$Relative.Growth.Rate+2*sqrt(fullSet2$variance)
	
	  
	comb <- merge(fullSet2, geneinfo, by.x = "gene", by.y = "Old.Gene.ID", all.x = TRUE)
	comb <- merge(comb,genomiclocations,by.x="current_version_ID",by.y="gene",all.x=TRUE)	
	
	comb$Confidence = -log(comb$variance)
	comb$Relative.Growth.Rate[!is.finite(comb$Confidence) ] = sample(1:1000/1000, nrow(comb[!is.finite(comb$Confidence), ]),replace=TRUE)
	comb$Confidence[!is.finite(comb$Confidence) ] = 0.1
	comb$Relative.Growth.Rate[comb$Confidence < 0.1 ] = sample(1:1000/1000, nrow(comb[comb$Confidence < 0.1, ]),replace=TRUE)
	comb$Confidence[comb$Confidence < 0.1 ] = 0.1

	comb<-mutate(comb, Confidence = ifelse(Confidence>10 , 10, Confidence))
	return(comb)
	}
gaussianMeanAndVariance2<- function(vals, variances){
		  df<-data.frame(value=vals,variance=variances)
		  df<-df[complete.cases(df),]
		 
		  vals=df$value 
		  variances=df$variance
		  if(length(vals)==1){
		  var<-variances[1]
		  mean<-vals[1]
		  }
		  else{
		 
		  precs=1/variances
		     
		  mean=sum(vals*precs)/sum(precs)
		 var1=(1/sum(precs))*(1/(length(vals)-1))*sum((vals-mean)**2/variances)
		var2=1/sum(precs) 
		if(is.na(var1)){var1<-0}
		var<-max(var1,var2)
		 }
		 return(data.frame(mean=mean,var=var))
	}

colsdf<-read.csv("defaultcolumns.csv")
expsdf<-read.csv("defaultexperiments.csv")
expsdf<-merge(experimentTable,expsdf,all.x=T)


shinyServer(function(input, output,clientData,session) {
ranges <- reactiveValues(x = NULL, y = NULL)
updateCheckboxGroupInput(session, "columns", choices=as.character(colsdf$column), selected = as.character(colsdf[colsdf$default=="yes",]$column))
#updateCheckboxGroupInput(session, "includedexperiments", )
choices<-c("abc")

    output$plot1 <- renderPlot({
    validate(
  need(!is.null(input$includedexperiments), "Loading, please wait..")
)
	if(input$defaultplot==TRUE){
	comb <- graphdata()
	if(input$HideUnfiltered){
	comb<-comb[comb$selected>0.5,]
	}
        toadd = NULL
        if (input$showErrorBars == TRUE) {
            toadd <- geom_errorbarh(data = otherdata, alpha=0.5, aes(xmin = xmin, xmax = xmax, y = thingie, x = xmin + xmax/2), 
                color = "gray")
        }
        
      
		withProgress(message="Drawing plot",value=0.65,{
		palettes=NULL 
		
if(input$xaxis=="--"){xaxis=NULL}else{xaxis=input$xaxis}
if(input$yaxis=="--"){yaxis=NULL}else{yaxis=input$yaxis}
if(input$fill=="--"){fill=NULL}else{fill=input$fill}
if(input$textcap=="--"){textcap=NULL}else{textcap=input$textcap}
if(input$color=="--"){color=NULL}else{color=input$color}

if(!is.null(textcap)){
	comb$textcap<-as.character(ifelse(as.character(comb$color)=="Unselected","",as.character(comb[,textcap])))
	}

if (is.null(yaxis)){
        p<-ggplot(comb, aes_string( x = xaxis, color = color,fill=fill, alpha = "selected"))  
			}
			else{
			    p<-ggplot(comb, aes_string(y = yaxis, x = xaxis, color = color,fill=fill, alpha = "selected"))  
			}
			p<-p +theme_gray(input$base_size)+palettes + 
            toadd + theme(legend.position = "bottom") +scale_alpha_continuous(limits=c(0,1),range=c(0,1)) + guides(alpha=FALSE) 
			
			if(input$plottype=="Point") {
			if(input$color=="color"){
			p=p+ geom_point(data=comb[comb$color=="Unselected",],size=input$pointsize,position=position_jitter(width=input$xjitter,height=input$yjitter))
			if(sum(comb$color!="Unselected")==1){
				p=p+ geom_point(data=comb[comb$color!="Unselected",],shape=1,size=20,color="black")
		
			}
			p=p+ geom_point(data=comb[comb$color!="Unselected",],size=input$pointsize,position=position_jitter(width=input$xjitter,height=input$yjitter))
			
			}
			else{
p=p+ geom_point(size=input$pointsize,position=position_jitter(width=input$xjitter,height=input$yjitter))
}
}

	if(input$plottype=="Violin") {
p=p+ geom_violin(color="blue") + stat_summary(fun.y = "median", colour = "red", size = 2, geom = "point")
}  

	if(input$plottype=="Density") {
p=p+ geom_density()
}

if(input$plottype=="Bar (mean)") {
p=p+  stat_summary(fun.y = "mean",  size = 2, geom = "bar")+   stat_summary(fun.data = mean_cl_normal, geom = "errorbar",position=position_dodge(), fun.args = list(mult = 1)) + expand_limits(y=0)
}
if(input$plottype=="Bar (identity)") {
p=p+  geom_bar(position="dodge")
}



if(input$plottype=="Box") {
p=p+  geom_boxplot()
}

if(input$xtrans!="identity") {
p=p+ scale_x_continuous(trans=input$xtrans)
}
if(input$ytrans!="identity") {
p=p+ scale_y_continuous(trans=input$ytrans)
}
if(!is.null(ranges$x))
{
p=p+ coord_cartesian(xlim = ranges$x, ylim = ranges$y)
}
else{
if(input$enablexlimit&input$enableylimit) {
p=p+coord_cartesian(xlim = c(input$xmin,   input$xmax),ylim = c(input$ymin,   input$ymax))
}
else{
			if(input$enablexlimit) {
p=p+coord_cartesian(xlim = c(input$xmin,   input$xmax))
}
if(input$enableylimit) {
p=p+coord_cartesian(ylim = c(input$ymin,   input$ymax))
}
}
}


if(input$rotatexlabels) {
p=p+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

if (input$color == "color") {
		
p=p+ scale_color_manual(values = phenolevelscolor)+labs(color="Growth phenotype")		
}

if (input$fill == "color") { 
		
p=p+ scale_fill_manual(values = phenolevelscolor)+labs(fill="Growth phenotype")		
}
if (input$xaxis == "Relative.Growth.Rate") {
		  
p=p+labs(x="Relative growth rate")		
} 
if(!is.null(textcap)){
p=p+geom_text(aes(label=textcap),hjust = 0 ) 
}


p
})
	} 
	else{
	 updateCheckboxInput(session, "HideMultiples",value = FALSE)
	 updateCheckboxInput(session, "mergemultipleobs",value = TRUE)
		con <- textConnection(input$customgraphinput)
		datastore <- read.table(con,header=FALSE)
		
		
		if(input$customtype=="number"){
	
		newdat<-merge(datastore,combSource(),by.x="V1",by.y=input$customformat)
		newdat$color=newdat$phenotype
		p<-ggplot(newdat,aes(x=as.numeric(as.character(V2)),y=color,color=color))+geom_point(alpha=0.5,position=position_jitter(width=0,height=0.4))+labs(x=input$customxlabel)
		#bladf<-data.frame(xa=1:10,ya=1:10)
		#ggplot(bladf,aes(x=xa,y=y))+geom_point()
		}
		if(input$customtype=="order"){
		datastore$V1=factor(as.character(datastore$V1),levels=as.character(datastore$V1))
		newdat<-merge(datastore,combSource(),by.x="V1",by.y=input$customformat,all.x=!(input$customhidemissing))
		if(length(newdat$V2)==0){newdat$V2=0;}
		
		p<-ggplot(newdat,aes(x=V1,y=Relative.Growth.Rate,ymin=lower,ymax=upper,fill=V2))+geom_bar(stat="identity")+geom_errorbar()+coord_cartesian( ylim = c(0, 1.5)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(x=input$customxlabel)
		}
		p+theme_bw()
	}
	})
    
	
	 graphdata <- reactive({
	
	 comb <- combSource()
	 
	 if(input$defaultplot==TRUE){
	  if (input$Pvalues == TRUE) {
            
            
            comb$color = comb$phenotype
            
        } else {
            comb$color = comb$barseq.phenotype
        }
		
		if (input$colorSelected == TRUE) {
            
            if(max(comb$selected>0.5)){
            comb$color = ifelse(comb$selected<0.5,"Unselected",as.character(comb$color))
			}
            
        }
		if (input$Pvalues == TRUE) {
			comb$color=factor(comb$color,levels=phenolevels)
		}
		
		 
	 }
	 if(input$xaxis=="bla"){currentx="Relative.Growth.Rate"}else{currentx=input$xaxis}
	  if(input$yaxis=="bla"){currenty="Confidence"}else{currenty=input$yaxis}
	  if(input$color=="bla"){currentcolor="color"}else{currentcolor=input$color}
	  if(input$fill=="bla"){currentfill="color"}else{currentfill=input$fill}
	 currenttextcap=input$textcap
	  curcolnames<-c("--",colnames(comb))
	  if(1  ){
	 updateSelectInput(session,"xaxis",choices=curcolnames,selected=currentx)
	   updateSelectInput(session,"yaxis",choices=curcolnames,selected=currenty)
	   updateSelectInput(session,"color",choices=curcolnames,selected=currentcolor)
	   updateSelectInput(session,"fill",choices=curcolnames,selected=currentfill)
	   updateSelectInput(session,"textcap",choices=curcolnames,selected=currenttextcap)
	   choices<<-curcolnames
	   }
	 return(comb)
	  
	 })
    # output$genetable <- renderText({ nearPoints(comb, input$plot1_click)
    
    
    # })
    output$call <- renderPrint({
    validate(
  need(storeRow2()$upper>0," ")
)
    upper<-storeRow2()$upper
    lower<-storeRow2()$lower
    if(upper<1&lower>0.1){
    cat("From this experiment, this mutant appears to be slow-growing. You might expect to see barcodes slowly decreasing in abundance over time in the plot below.")
    }
     if(upper>1&lower>0.1){
    cat("From this experiment, this mutant may be dispensable. You might expect to see barcodes increasing over time, or remaining at a relatively constant level, in the plot below.")
    }
    
    if(upper<1&lower<0.1){
    cat("From this experiment, this mutant may be essential. You might expect to see barcode abundance falling rapidly in the plot below.")
    }
     if(upper>1&lower<0.1){
    cat("From this experiment, it was not possible to call a phenotype for this mutant. There may be inconistent results, or the mutant may establish at such a low level that there a very few readings.")
    }
    
    #cat("hi")
    })
	oldsel<-c(0)
    storeRow <- reactive({
    
    validate(
  need(!is.null(input$includedexperiments), "Loading, please wait..")
)
if( sum(combSource()$selected>0.5)==1){
return(combSource()[combSource()$selected>0.5,]);
}
        if(input$mainPanel=="plot")
	{
if( sum(combSource()$selected>0.5)>0){
	
	return(nearPoints(combSource()[combSource()$selected>0.5,], input$plot1_click)[1, ])
	}
	else{
	return(nearPoints(combSource(), input$plot1_click)[1, ])
	}
	
		}
		
		else{
		
		newsel<-input$maintable_rows_selected
		if (length(newsel)>0){
		cs<-combSource()
		start<-cs[cs$selected>0.5,]
		return(start[newsel,])
		} 
		else{
		bla<-data.frame(gene=NA)
		return(bla)
		
		}
		
		}
    }) 
	
	storeRow2 <- reactive({
	rows<-singlecomb()
        rows[rows$gene==storeRow()$gene & rows$experiment==input$moredetailsselect,]
    })
    
    combSource <- reactive({
	withProgress(message="Loading secondary data",value=0.5,{
	
        newcomb <- initcomb()
        
		#tab <- as.data.frame(table(newcomb$gene))
       ###START THE BLOCK CONTAINING THE PROCESS OF WORKING OUT WHAT IS SELECTED
        newcomb$selected = input$pointopacity/100
        if (input$genelist != "") {
        
            glist = unlist(strsplit(input$genelist, "\n"))
            glist<- sapply(glist,function (x) gsub("^\\s+|\\s+$", "", x))
            newcomb$selected = newcomb$selected * (as.numeric(newcomb$gene %in% glist | newcomb$current_version_ID %in% glist | newcomb$PfID %in% glist))
            
        } 
		
		
		
		
		if (input$search != "") {
            
            
            newcomb$selected = newcomb$selected * (as.numeric(grepl(input$search, newcomb$gene_product, ignore.case = TRUE) | 
                grepl(input$search, newcomb$gene_name, ignore.case = TRUE)|grepl(input$search, newcomb$gene, ignore.case = TRUE)|grepl(input$search, newcomb$current_version_ID, ignore.case = TRUE)|grepl(input$search, newcomb$cloneid, ignore.case = TRUE)))
            
            
        }
		if (length(input$enrichmenttable_rows_selected) >0) {
		goterms<-enrichmentReactive()$table[input$enrichmenttable_rows_selected,1]
		if(length(goterms>0)){
		godata<-enrichmentReactive()$godata
		g<-unlist(genesInTerm(godata,goterms))
		newcomb$selected = newcomb$selected * (as.numeric(newcomb$current_version_ID %in% g | newcomb$PfID %in% g))
		}
		}
        if (input$filterconf ==TRUE) {
         newcomb$selected = newcomb$selected * (as.numeric(newcomb$Confidence>=input$confidence[1]&newcomb$Confidence<=input$confidence[2]))
        }
			if (input$filterfitness ==TRUE) {
            
            if (input$filterfitnesstype=="bestguess"){
            newcomb$selected = newcomb$selected * (as.numeric(newcomb$Relative.Growth.Rate>input$fitness[1]&newcomb$Relative.Growth.Rate<input$fitness[2]))
			}
			else{
			lowbar=input$fitness[1]
			highbar=input$fitness[2]
			if(lowbar==0){lowbar=-Inf}
			if(highbar==1.5){highbar=Inf}
			newcomb$selected = newcomb$selected * (as.numeric(!is.na(newcomb$lower)&!is.na(newcomb$upper)&newcomb$lower>lowbar&newcomb$upper<highbar))
			}
            #&newcomb$upper>input$fitnessupper[1]&newcomb$upper<input$fitnessupper[2]&newcomb$lower>input$fitnesslower[1]&newcomb$lower<input$fitnesslower[2]
            
        }
		exps<-input$experiments_rows_selected
		if (length(exps)>0 & input$mergemultipleobs==FALSE) {
		valid=experimentTable[exps,"Experiment"]
		newcomb$selected = newcomb$selected * (newcomb$experiment %in% valid)
		
		}
        newcomb[newcomb$selected<0.1,'selected']= 0.25
        ##END THE BLOCK WORKING OUT WHAT IS SELECTED
        
        #newcomb<-newcomb[newcomb$phenotype!=phenolevels[1],]
        addPhenotypes(newcomb)
		})
    })
    
    
    library(arrayhelpers)
    output$plot2 <- renderPlot({
        row <- storeRow2()
		
        fileno = which(files == row$file)
        arrayforfile = arrays[[fileno]]
        id = which(arrayforfile$genes == row$gene)
        
        ratios <- array2df(arrayforfile$ratios[, , id], label.x = "ratio", levels = list(mouse = NA, day = NA))
        ratiosvar <- array2df(arrayforfile$ratios[, , id], label.x = "ratiovar", levels = list(mouse = NA, day = NA))
        merge <- cbind(ratios, ratiosvar)
        merge$sd = sqrt(merge$ratiovar)
        merge$ratiomin = merge$ratio - merge$sd * 2
        merge$ratiomax = merge$ratio + merge$sd * 2
        ggplot(merge, aes(x = day + 3, y = ratio, color = as.factor(mouse))) + geom_point() + geom_line() + ggtitle(row$gene) + 
            expand_limits(y = 0) + labs(x = "Day", color = "Mouse")
        
    })
    
    output$plot3 <- renderPlot({
        row <- storeRow2()
        validate(
  need(nchar(row$file)>0," ")
)
        fileno = which(files == row$file)
        arrayforfile = arrays[[fileno]]
        id = which(arrayforfile$genes == row$gene)
        
        abs <- array2df(arrayforfile$absfitness[, , id], label.x = "abs", levels = list(mouse = NA, day = NA))
        absvar <- array2df(arrayforfile$absfitnessvar[, , id], label.x = "absvar", levels = list(mouse = NA, day = NA))
        merge <- cbind(abs, absvar)
        merge$sd = sqrt(merge$absvar)   
        merge$absomin = merge$abs - merge$sd
        merge$absmax = merge$abs + merge$sd
        ggplot(merge, aes(x = day + 4, y = abs, color = as.factor(mouse))) + geom_bar(stat = "identity", position = "dodge") + 
            ggtitle(row$gene) + expand_limits(y = 0) + labs(x = "Day", color = "Mouse",y="Normalised relative growth rate")+coord_cartesian(ylim=c(0,2))+ geom_hline(yintercept=1, color="black",alpha=0.15,size=2)
        
    })
    
    output$plot4 <- renderPlot({
        row <- storeRow2()
        fileno = which(files == row$file)
        arrayforfile = arrays[[fileno]]
        id = which(arrayforfile$genes == row$gene)
        
        abs <- array2df(arrayforfile$absfitness[, , id], label.x = "abs", levels = list(mouse = NA, day = NA))
        absvar <- array2df(arrayforfile$absfitnessvar[, , id], label.x = "absvar", levels = list(mouse = NA, day = NA))
        merge <- cbind(abs, absvar)
        merge$inversevar = 1/merge$absvar
        ggplot(merge, aes(x = day + 4, y = inversevar, fill = as.factor(mouse))) + geom_bar(position = "dodge", stat = "identity") + 
            ggtitle(cat("Fitness: ", row$gene)) + scale_y_continuous(expand = c(0, 0)) + labs(x = "Day", fill = "Mouse",y="Weight (precision)")
         
    })
	 output$plotcontrols <- renderPlot({
        row <- storeRow2()
        fileno = which(files == row$file)
        arrayforfile = arrays[[fileno]]
       cnames<-arrayforfile$controlnames
        
        control <- array2df(arrayforfile$controlarray[, , ], label.x = "control", levels = list(mouse = NA, day = NA,gene=NA))
        controlvar <- array2df(arrayforfile$controlvararray[, , ], label.x = "controlvar", levels = list(mouse = NA, day = NA,gene=NA))
        merge <- cbind(control, controlvar)
        merge$inversevar = 1/merge$controlvar
		merge$gene<-cnames[merge$gene]
        ggplot(merge, aes(x = day + 4, y = inversevar, fill = as.factor(mouse))) + geom_bar(position = "dodge", stat = "identity") + 
            ggtitle(cat("Control genes: ", row$gene)) + scale_y_continuous(expand = c(0, 0)) + labs(x = "Day", fill = "Mouse")+facet_grid(.~gene)
         
    })
    output$genetable <- renderUI({
        # Because it's a ggplot2, we don't need to supply xvar or yvar; if this were a base graphics plot, we'd need those.
		toattempt<-c("Relative.Growth.Rate", "lower", "upper",  "experiment","call","timesAnalysed")
		rower<-storeRow() 
		present<-which( colnames(rower)%in% toattempt)
        bla <- storeRow()[, present]
		
		bla[,sapply(bla, is.numeric)] <-round(bla[,sapply(bla, is.numeric)],2)
        df = data.frame(b = t(bla))
		
        #print(df)
    
       string<-paste("<table><tr><td class='lefttd'>Relative growth rate:</td><td>",round(rower$Relative.Growth.Rate,2),"</td>","</tr><tr><td class='lefttd'>95% CI:</td><td>",round(rower$lower,2)," – ",round(rower$upper,2),"</td>",sep="",collapse="");
          string<-paste(string,"<tr><td class='lefttd'>Phenotype call:</td><td style='color:",phenolevelscolor[rower$phenotype],"'>",phenolevelslong[rower$phenotype],"</td>","</tr>",sep="",collapse="");
        
        if ("timesAnalysed" %in% colnames(rower)){
        string<-paste(string,"<tr><td class='lefttd'>Times analysed:</td><td>",rower$timesAnalysed,"</td>","</tr>",sep="",collapse="");
        }
        if ("experiment" %in% colnames(rower)){
        string<-paste(string,"<tr><td class='lefttd'>Experiment:</td><td>",rower$experiment,"</td>","</tr>",sep="",collapse="");
        }
        
        cloneidnum=gsub("[^0-9]","",rower$cloneid)

 		if(!is.na(cloneidnum) & !is.na(rower$cloneid) & cloneidnum != "" ){
 		 string<-paste(string,"<tr><td class='lefttd'>Clones transfected:</td><td><a href='http://plasmogem.sanger.ac.uk/designs/final_vector/", cloneidnum,"'>",rower$cloneid,"</a></td>","</tr>",sep="",collapse="");
     
 		}
        string<-paste(string,"</table>",sep="",collapse="");
        
     HTML(string)
      
        
    })
    
        output$genefitness <- renderPlot({
        rower<-storeRow() 
        df=data.frame(x=seq(-1,1.5,length=500))
        df$y=dnorm(df$x,rower$Relative.Growth.Rate,sqrt(rower$variance))
        ggplot(df,aes(x=x,y=y))+geom_line(color=phenolevelscolor[rower$phenotype])+ coord_cartesian(xlim = c(0, 1.25))+labs(x="Relative Growth Rate",y="Likelihood")+theme(axis.ticks.y=element_blank())+scale_x_continuous( breaks = c(0,0.5,1))+ scale_y_continuous(breaks=NULL)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme_bw()
        
    })
	
	round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))

  df[,nums] <- round(df[,nums], digits = digits)

  (df)
}
	output$genetable2 <- renderTable({
        # Because it's a ggplot2, we don't need to supply xvar or yvar; if this were a base graphics plot, we'd need those.
		toattempt<-c("Relative.Growth.Rate", "lower", "upper","experiment","phenotype")
		single<-singlecomb()
		present<-which(colnames(single)%in%toattempt);
        gene <- storeRow()[, "gene"]
        if("phenotype"%in%colnames(single)){
		 bla <- single[single$gene==gene, c("experiment","Relative.Growth.Rate", "lower", "upper","phenotype")]
		 }
		 else{
		 bla <- single[single$gene==gene, c("experiment","Relative.Growth.Rate", "lower", "upper")]
		 }
		 bla<-round_df(bla, digits=2)
		 #experiment<-storeRow()[,"experiment"]
	
		 if("experiment"%in% colnames(storeRow())){
		 updateSelectInput(session, "moredetailsselect", choices=bla$experiment,selected=storeRow()[, "experiment"])
		 }
		 else{
		 
	
		 updateSelectInput(session, "moredetailsselect", choices=bla$experiment)
		 }
       print(bla)
        
    },include.rownames=FALSE)
    
    output$geneid <- renderPrint({
	if(is.na(as.character(storeRow()$gene))){
	cat("")
	}
	else{
	if(input$newGeneID){
	cat(as.character(storeRow()$current_version_ID))
	}
	else{
	cat(as.character(storeRow()$gene))
	}
        }
    })
    
    output$genename <- renderPrint({
        cat(as.character(storeRow()$gene_name))
    })
    
    output$geneproduct <- renderPrint({
        cat(as.character(storeRow()$gene_product))
    })
    output$searchResult <- renderPrint({
	cat ("Displaying a total of ")
	cat(length(combSource()$selected))
	cat (" values, of which ")
        cat(sum(combSource()$selected>0.5))
	cat (" pass the filtering conditions.")	
    })
    output$searchResult2 <- renderPrint({
    validate(
  need(!is.null(input$includedexperiments), "")
)
    if(sum(combSource()$selected>0.5)==0){cat("No phenotyped genes satisfy your filters.")}

    })
	
	  output$downloadwarning <- renderPrint({
	  selected=sum(combSource()$selected>0.5)
	  total=length(combSource()$selected)
	cat ("A CSV file with ")
        cat(selected)
	cat (" rows will be produced.")
if(selected<total){
cat (" Please note that this is not the full dataset, as you have entered filtering conditions for the data.")
}	
    })
    
    output$enrichmenttable <- renderDataTable({
	df<-enrichmentReactive()$table
	colnames(df)=c("GO ID","Term","Total (here)","Selected","Expected","p")
	df
	}
	)
      output$columns <- renderDataTable({

	}
	)
    
    goDBReactive <- reactive({
        if (input$gospecies == "pb") {
            gomap <- read.csv("./otherdata/geneid2gopb.csv", stringsAsFactors = FALSE)
            
        } else {
            gomap <- read.csv("./otherdata/geneid2gopf.csv", stringsAsFactors = FALSE)
        }
        go <- aggregate(GO ~ ID, data = gomap, c)
        
        godb <- setNames(as.list(go$GO), go$ID)
        godb
    }) 
	 
	goDBReactive2 <- reactive({
        if (input$gospecies == "pb") {
            gomap <- read.csv("./otherdata/geneid2gopb.csv", stringsAsFactors = FALSE)
            
        } else {
            gomap <- read.csv("./otherdata/geneid2gopf.csv", stringsAsFactors = FALSE)
        }
       return(gomap)
    })
    
    enrichmentReactive <- reactive({
        withProgress(message = "Loading GO data", value = 0.15, {
		comb <- isolate(combSource()) 
		input$goGO
		if(length(isolate(input$enrichmenttable_rows_selected))==0){
            
            godb <- isolate(goDBReactive())
            if (isolate(input$gospecies) == "pb") {
                allgenes = setNames(rep(1, length(comb$current_version_ID)), comb$current_version_ID)
                myInterestingGenes <-isolate( brushedPoints(comb, input$plot1_brush)$current_version_ID)
                geneList <- factor(as.integer(comb$current_version_ID %in% myInterestingGenes))
                names(geneList) <- comb$current_version_ID
                
                
            } else {
                allgenes = setNames(rep(1, length(comb$PfID)), comb$PfID)
                myInterestingGenes <- isolate(brushedPoints(comb, input$plot1_brush)$PfID)
                geneList <- factor(as.integer(comb$PfID %in% myInterestingGenes))
                names(geneList) <- comb$PfID 
                
            }
            validate(
  need(length(levels(geneList))==2, "Please drag a box in the plot to select some genes")
)
            
            
            isolate({
            GOdata <- new("topGOdata", ontology = input$ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = godb)
            setProgress(message = "Calculating statistics", value = 0.5)
            resultant <- runTest(GOdata, algorithm = input$goalgo, statistic = input$gostat)
            setProgress(message = "Calculating statistics", value = 0.65)
            allRes <- GenTable(GOdata, res = resultant, orderBy = "res", ranksOf = "res", topNodes = min(200,length(resultant@score)))
        })}
		})
        
        list(table=allRes,godata=GOdata);
    })
    
     
     output$includedrender<- renderUI({
     checkboxGroupInput("includedexperiments", "Included experiments",choices=as.character(expsdf$Experiment), selected = as.character(expsdf[!is.na(expsdf$default)&expsdf$default=="yes",]$Experiment))
     }
     );
     outputOptions(output, "includedrender", suspendWhenHidden = FALSE)
     outputOptions(output, "includedrender", priority = 10)
    output$downloadData <- downloadHandler(filename = function() {
    if(input$downloadfilename==""){
        return(paste("Barseq", format(Sys.Date(), format = "%Y%m%d"), ".csv", sep = ""))
        }
       else {
        return(paste(input$downloadfilename, ".csv", sep = ""))
        }
        
    }, content = function(file) {
	cs<-combSource()
	cs<-cs[cs$selected>0.5,]
	cs<-cs[,input$columns[input$columns %in% colnames(cs)]]
        write.csv(cs, file,row.names=F) 
    })
	
	
	singlecomb<-reactive({ 
	if (is.null(input$includedexperiments) || is.na(input$includedexperiments)){
  return()
}
withProgress(message="Loading initial data",value=0.25,{
if(input$usermgmdb==TRUE){
	fullSet2 <- read.csv("./otherdata/rmgmdbinfakeformat.csv")
	fullSet2$variance=sample((1:1000/20000),length(fullSet2$variance))
	}else{
	if(input$HideTags==FALSE){
  fullSet2<-fullSet
  }
  else{
   
   fullSet2<-fullSet[!(grepl("tag", fullSet$gene, ignore.case = TRUE)),]
   
  }
	}

	fullSet2<-fullSet2[fullSet2$experiment %in% input$includedexperiments,] 
	

	
	fullSet2<-addExtraData(fullSet2)
	})
	#cat("hello!", file = stderr())
	#fullSet2
	#fullSet2
	fullSet2<-merge(fullSet2,barseqtemp,by=c("experiment"),all.x=TRUE)
	fullSet2<-merge(fullSet2,cloneIDs,by=c("gene","barseq"),all.x=TRUE)
fullSet2<-merge(fullSet2,homology,by=c("cloneid"),all.x=TRUE)

	
	
})

 addPhenotypes<-function(newcomb){
#NOW CALCULATE PHENOTYPES USING P VALUES
        newcomb$z1 <- (1 - newcomb$Relative.Growth.Rate)/sqrt(newcomb$variance)
        newcomb$p1 <- 2 * pnorm(-abs(newcomb$z1))
        newcomb$f1 = p.adjust(newcomb$p1, method = input$padjmethod)
        newcomb$call1 = FALSE
        newcomb[newcomb$f1 < input$pvalue, ]$call1 = TRUE
        
        newcomb$z0 <- (0.1 - newcomb$Relative.Growth.Rate)/sqrt(newcomb$variance)
        newcomb$p0 <- pnorm(-abs(newcomb$z0))
        newcomb$f0 = p.adjust(newcomb$p0, method = input$padjmethod)
        newcomb$call0 = FALSE
        newcomb[newcomb$f0 < input$pvalue & newcomb$Relative.Growth.Rate > 0.1, ]$call0 = TRUE
        newcomb$phenotype = "None"
        newcomb$phenotype[newcomb$call0 & newcomb$call1&newcomb$Relative.Growth.Rate>1 ] <- phenolevels[5]
		newcomb$phenotype[newcomb$call0 & newcomb$call1 &newcomb$Relative.Growth.Rate<1] <- phenolevels[3]
        newcomb$phenotype[newcomb$call0 & !newcomb$call1 ] <- phenolevels[4]
        newcomb$phenotype[!newcomb$call0 & newcomb$call1 ] <- phenolevels[2]
        newcomb$phenotype[!(newcomb$call0 | newcomb$call1) ] <- phenolevels[1]
		newcomb$phenotype=factor(as.character(newcomb$phenotype),levels=phenolevels)
      
        newcomb$id=1:length(newcomb$phenotype)
        return(newcomb);
}
multicomb<-reactive({ 


	fullSet2<-singlecomb()
	if(is.null(fullSet2)){return(NULL);}
	fullSet3<- fullSet2 %>% group_by(gene) %>% do(gaussianMeanAndVariance2(.$Relative.Growth.Rate,.$variance)) %>% transmute(Relative.Growth.Rate=mean, variance=var)
	fullSet4<-fullSet2  %>%  group_by(gene)  %>%  summarise(cloneid=paste(unique(cloneid),sep=",",collapse=","),timesAnalysed=length(Relative.Growth.Rate),normd6toinputA=mean(normd6toinputA,na.rm=T),normd6toinputB=mean(normd6toinputB,na.rm=T),normd6toinputC=mean(normd6toinputC,na.rm=T))

	fullSet2<-merge(fullSet3,fullSet4,by=c("gene")) 
	fullSet2<-addExtraData(fullSet2)
	fullSet2
	})
initcomb<-reactive({ 
if(input$mergemultipleobs==TRUE){
	multicomb()
}
else{
singlecomb()}
})


 observeEvent(input$egTagKO, {
 print("bla")
   fileName <- './otherdata/tagvskos.txt'
chars<-readChar(fileName, file.info(fileName)$size)
updateTextInput(session,inputId ="customgraphinput",value=chars)
updateSelectInput(session,"customtype",selected="order")
updateSelectInput(session,"customformat",selected="gene")
updateTextInput(session,"customxlabel",value="Gene")
updateCheckboxInput(session,"HideTags",value=FALSE)
  })
  
   observeEvent(input$egAP2O, {
   fileName <- './otherdata/ap2oagged.txt'
chars<-readChar(fileName, file.info(fileName)$size)
updateTextInput(session,inputId ="customgraphinput",value=chars)
updateSelectInput(session,"customtype",selected="number")
updateSelectInput(session,"customformat",selected="gene")
updateTextInput(session,"customxlabel",value="log2 fold change in AP2-O mutant")
  })
  
  observeEvent(input$egRibosome, {
 print("bla")
   fileName <- './otherdata/ribosomevstransport.txt'
chars<-readChar(fileName, file.info(fileName)$size)
updateTextInput(session,inputId ="customgraphinput",value=chars)
updateSelectInput(session,"customtype",selected="order")
updateSelectInput(session,"customformat",selected="current_version_ID")
updateTextInput(session,"customxlabel",value="Gene")
updateCheckboxInput(session,"customhidemissing",value=TRUE)

  }) 
output$experiments <- DT::renderDataTable(experimentTable)

#output$maintable <- DT::renderDataTable({iris}
#)
output$maintable <- DT::renderDataTable({
cs<-combSource()
if(input$newGeneID){
	start<-cs[cs$selected>0.5,c("current_version_ID","gene_name","gene_product","Relative.Growth.Rate","lower","upper","phenotype")]
	}
	else{
	start<-cs[cs$selected>0.5,c("gene","gene_name","gene_product","Relative.Growth.Rate","lower","upper","phenotype")]
	}

start$upper=round(start$upper,2)
start$lower=round(start$lower,2)
start$Relative.Growth.Rate=round(start$Relative.Growth.Rate,2)
start
},selection="single",rownames=FALSE,colnames=c("Gene ID","Name","Product","Growth Rate","lower","upper","Phenotype")
)

output$debug <- renderPrint({


})
    
      

 output$selectedgeneinfo <- renderPlot({
 cs<-combSource()
 cs$searched="Not in selection"
cs[cs$selected>0.5,]$searched="In selection"
 ggplot(cs,aes(x=phenotype,fill=phenotype))+geom_bar()+facet_wrap(~searched,scales="free")+ theme(legend.position = "bottom")+scale_fill_manual(values = phenolevelscolor)	+guides(fill=guide_legend(nrow=2,byrow=TRUE))+ theme( axis.text.x = element_blank())
 })
output$knitDoc <- renderUI({
    input$cceval
    return(isolate(HTML(knit2html(text = input$ccrmd, fragment.only = TRUE, quiet = TRUE))))
  })  
  output$report = downloadHandler(
    filename = 'myreport.pdf',
    
    content = function(file) {
     out <- tempfile(fileext='.md')
    knit(text = input$ccrmd, output=out)
     out <- pandoc(input=out, format='latex')
      file.rename(out, file) # move pdf to file for downloading
    },
    
    contentType = 'application/pdf'
  )
 filesreac<-reactive({
input$ccrefresh
 return(data.frame(file=list.files("./customanalyses/", full.names = FALSE)))
 })
 
 output$ccfiles <- DT::renderDataTable({
 filesreac()
 },selection="single")
 observeEvent(input$ccsave, {
    write(input$ccrmd, file = paste("./customanalyses/", input$ccfilename, sep=""))
  })
  
  
   observeEvent(input$ccload, {
   filesel<-input$ccfiles_rows_selected
   files<-filesreac()
   fileName <- paste("./customanalyses/", filesreac()[filesel,"file"], sep="")
chars<-readChar(fileName, file.info(fileName)$size)
updateAceEditor(session, "ccrmd", chars,   mode="r")
  })
  output$help <- renderText({
  
 fileName <- "./help.html"
chars<-readChar(fileName, file.info(fileName)$size)
print(chars)
#print("woohoo")
 })
#
 observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)

    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

 shinyURL.server()
})

