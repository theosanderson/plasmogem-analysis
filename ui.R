library(shiny)
library(shinyBS)
library(DT)
library(shinyAce)
admin=T
library("shinyURL")
mycss <- "
#plot-container {
  position: relative;
}
#loading-spinner {
  position: absolute;
  left: 50%;
  top: 50%;
  z-index: -1;
  margin-top: -50px;  /* half of the spinner's height */
  margin-left: -50px; /* half of the spinner's width */
}
#plot1.recalculating {
  z-index: -2;
  }
  body {
    overflow-y: scroll;
}

#debug{display:none}
#maintable{margin-bottom:20px;}
div#geneside{margin-bottom:30px;}

#searchResult{color:gray;font-size:80%;}
#searchResult2{color:red;font-weight:bold;}
.lefttd{text-align:right;padding-right:10px; font-weight:bold}
div#genetable{margin-top:10px; margin-bottom:10px;}

"
shinyUI(fluidPage(titlePanel("", windowTitle = "PlasmoGEM data explorer"),
  tags$head(tags$style("#genetable{color: black;
                       font-size: 12px;word-wrap:normal
                       }"
                         )),  tags$head(tags$style(HTML(mycss))),HTML('<div style="margin-top:15px;margin-bottom:15px;margin-left:5px; margin-right:5px;" id="header">
        <!-- this is the displayed page header, not the part between <head></head> tags -->
<a href="javascript:location.reload()"><img style="float:left" id="plasmogemLogo" src="http://theo.io/shiny/plasmogembeta.png" ,="" alt="PlasmoGEM logo" title="PlasmoGEM beta release, data liable to change"></a>
<img style="float:right" id="sangerLog" src="http://plasmogem.sanger.ac.uk/static/images/sanger-logo.png" ,="" alt="Sanger logo">
<div style="clear:both"></div>

      </div>'),HTML(" <script type='text/javascript'>
        $(function()
        {
            // Prevent accidental navigation away
            $(':input').bind(
                'change', function() { setConfirmUnload(true); });
            $('.noprompt-required').click(
                function() { setConfirmUnload(false); });

            function setConfirmUnload(on)
            {
                window.onbeforeunload = on ? unloadMessage : null;
            }
            function unloadMessage()
            {
                return ('You have entered new data on this page. ' +
                        'If you navigate away from this page without ' +
                        'first saving your data, the changes will be lost.');
            }

            window.onerror = UnspecifiedErrorHandler;
            function UnspecifiedErrorHandler()
            {
                return true;
            }

        }); 

    </script>"),

  fluidRow(column(9,bsModal("modalURL", "More information", "modalURLShow", size = "small",shinyURL.ui()),
  bsModal("modalColumns", "Customise columns", "modalColumnsShow", size = "large",
  checkboxGroupInput("columns", "Included columns", "")
  
  ),
  bsModal("modalExperiments", "Included experiments", "includedExperimentsShow", size = "large",
  htmlOutput("includedrender")
  
  ),
   bsModal("modalDetails", "More information", "modalDetailsToggle", size = "large",
   h3("Individual experiment results:"),
   tableOutput(outputId="genetable2"),HTML("<hr>"),
selectInput("moredetailsselect", "Display experiment:",choices=c(" ")),


h3("Plot of barcode ratios in the population"),textOutput("call"),
plotOutput("plot2",height="300px"),textOutput("inforatio"),h3("Normalised growth rate at each timepoint/mouse"),
plotOutput("plot3",height="300px"),
HTML("The growth rate values shown in the plot above are all estimates for the true growth rate of the mutant. But some of them are more accurate than others. The plot below shows the algorithm's estimate for the accuracy of each value. The most accurate values are weighted higher in calculating the final estimate for the mutant's growth rate."),
h3("Estimated accuracy of fitnesses at each timepoint/mouse"),plotOutput("plot4",height="300px")
#,plotOutput("plotcontrols",height="300px")
),

     
  tabsetPanel(type = "tabs",id="mainPanel", tabPanel("Plot",value="plot",
  div(id = "plot-container",
        tags$img(src = "http://theo.io/shiny/ajax-loader.gif",
                 id = "loading-spinner"),
         plotOutput("plot1", height = "400px", 
                             click = "plot1_click",
        dblclick = "plot1_dblclick",
        brush = brushOpts(
          id = "plot1_brush",
          resetOnNew = TRUE
        ))
    )
                 ),tabPanel("Table",value="table",
                  DT::dataTableOutput("maintable"
                            )))),
           
           column(3,
                  tabsetPanel(type = "tabs", id="topright",
                              tabPanel("Gene info",value="geneinfo",
                                       h4(textOutput("geneid")),div(id="geneside",conditionalPanel(condition = "output.geneid == ''",div("Click on a gene in the fitness plot for further information",style="opacity:50%;font-style:italic;")),conditionalPanel(condition = "output.geneid != ''",h5(textOutput("genename")),h5(textOutput("geneproduct")),htmlOutput(outputId="genetable"),plotOutput(width="100%", height="100px","genefitness"),
									   
									   actionButton("modalDetailsToggle", "More diagnostic info"))),
									   #actionLink("modalURLShow", "Copy a link to the current analysis"),
									   textOutput("searchResult"), textOutput("searchResult2")),
							tabPanel("Selected genes",
                                      plotOutput("selectedgeneinfo"))
									   
                              
                  ))),
  
  fluidRow(column(12,
                  tabsetPanel(type = "tabs", id="tooltab",
                              tabPanel("Text search",
                                       
                                       
                                       textInput("search", "Search (e.g. 'ribosomal' or 'ookinete' or 'CRT' or 'PBANKA_010490'):", "")),
                              tabPanel("Gene ID search",
                                       HTML('Gene list:<br><textarea id="genelist" rows="3" cols="40"></textarea>')
                                      
                              ),
							   tabPanel("Experiments",actionButton("includedExperimentsShow", "Included experiments"),  checkboxInput(inputId = "mergemultipleobs",  label = strong("Merge multiple observations"),value=TRUE  ),   conditionalPanel(condition = "input.mergemultipleobs == false",DT::dataTableOutput('experiments'))),
							   tabPanel("Relative growth rate",checkboxInput(inputId = "filterfitness",  label = strong("Filter on relative growth rate"),value=FALSE  ),   conditionalPanel(condition = "input.filterfitness == true",selectInput("filterfitnesstype", "Type:",c("Best guess" = "bestguess", "Confidence interval" = "confint")),sliderInput("fitness","Relative Growth Rate",0,1.5,c(0,1.5))),
							   checkboxInput(inputId = "filterconf",  label = strong("Filter on confidence"),value=FALSE  ),   conditionalPanel(condition = "input.filterconf == true",sliderInput("confidence","Fitness",0,10,c(0,10)))
							   ),
                             
                              tabPanel("Gene Ontology enrichment",
                                        fluidRow(column(3,selectInput("gospecies", "Species:",c("P. berghei" = "pb", "P. falciparum" = "pf"))),
                                                 column(3,selectInput("ontology", "Ontology:",c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"))),column(3,checkboxInput("showAdvancedGO","Show advanced options",value=FALSE))), actionButton("goGO", "Calculate"),
                                        conditionalPanel(condition = "input.showAdvancedGO == true",
                                                         selectInput("goalgo", "Algorithm:",c( "weight01"= "weight01","Classic" = "classic","elim"="elim","weight"="weight" ,"lea"="lea","Parent-child"="parentchild"),selected="classic"),
                                                         selectInput("gostat", "Statistic:",c("Fisher" = "fisher", "KS"= "ks","T"= "t", "Global test" ="globaltest", "Sum"= "sum")
                                                                     
                                                         )),dataTableOutput("enrichmenttable")
                              )
                              ,
							  tabPanel("Advanced settings",
							   tabsetPanel(type = "pills",
							   
									  tabPanel("Visualisation",checkboxInput(inputId = "newGeneID",  label = strong("Use new gene IDs")),checkboxInput(inputId = "showErrorBars",  label = strong("Error guides")  ),checkboxInput(inputId = "colorSelected",  label = strong("Colour selection"),value=TRUE  ),  sliderInput("pointopacity","Point opacity",0,100,70),sliderInput("pointsize","Point size",0,2,1.3,step=0.1)),
									   tabPanel("Phenotypes",
                                       checkboxInput(inputId = "Pvalues",  label = strong("Calculate new phenotypes") ,value=TRUE ),
                                       selectInput("padjmethod", "Multiple comparison adjustment:",
                                                   c("FDR" = "fdr", "None" = "none", "bonferroni" = "bonferroni")),
                                       sliderInput("pvalue","p-value",0,1,0.05) 
                              )
                                      ,tabPanel("Custom graph",sliderInput("base_size","Scaling",1,30,15),selectInput("xaxis",choices=c("Relative.Growth.Rate"),label="X"),selectInput("yaxis",choices=c("Confidence"),label="Y"),selectInput("color",choices=c("color"),label="Color"),selectInput("fill",choices=c("--"),label="Fill"),selectInput("textcap",choices=c("--"),label="Text"),selectInput("xtrans",choices=c("identity","log10"),label="X trans"),selectInput("ytrans",choices=c("identity","log10"),label="Y trans"),checkboxInput("enablexlimit",label="X limits",value=TRUE),
					conditionalPanel(condition = "input.enablexlimit == true",numericInput("xmin",label="X min",value=0),numericInput("xmax",label="X max",value=1.3)),checkboxInput("enableylimit",label="Y limit"),
					conditionalPanel(condition = "input.enableylimit == true",numericInput("ymin",label="Y min",value=0),numericInput("ymax",label="Y max",value=10)),selectInput("plottype",label="Plot type",choices=c("Point","Bar (mean)","Bar (identity)","Violin","Density","Box")),checkboxInput(inputId = "rotatexlabels",  label = strong("Rotate X-axis labels")),checkboxInput(inputId = "HideUnfiltered",  label = strong("Show only items matching search")),numericInput("xjitter",label="X jitter",value=0),numericInput("yjitter",label="Y jitter",value=0),checkboxInput(inputId = "defaultplot",  label = strong("Default plot"),value=TRUE  ),conditionalPanel(condition = "input.defaultplot == false",
							  actionButton("egTagKO", "Example 1: tags vs KOs"),actionButton("egAP2O", "Example 2: AP2-O expression level"),actionButton("egRibosome", "Example 3: Ribosome vs. transporter"),
							  HTML('<br>Tabular input (first column = gene, optional second column = value/category ):<br><textarea id="customgraphinput" rows="3" cols="40"></textarea>'),
selectInput("customformat", "Gene format",
                                                   c("Pb old" = "gene", "Pb new" = "current_version_ID", "Pf new" = "PfID"))
,
selectInput("customtype", "Table type",
                                                   c("Gene and number" = "number", "Genes in order" = "order"))
,							  plotOutput("plotcustom"),
textInput("customxlabel","X-axis label"),
checkboxInput("customhidemissing","Hide missing values"))

                              ),
                              
                              tabPanel("Pre-processing",checkboxInput(inputId = "usermgmdb",  label = strong("Use RMGMDB instead"),value=FALSE  ),checkboxInput("addpfextra",label="Add extra data from Pf orthologues"),checkboxInput(inputId = "HideTags",  label = strong("Hide tags"),value=TRUE  )
									  ))
							  ),
							  
									   
                                       
									    
										 
                                     
                             
                              
                               tabPanel("Download CSV",  textOutput("downloadwarning"), textInput("downloadfilename","File name (optional)"),actionButton("modalColumnsShow", "Customise columns"),downloadButton('downloadData', 'Download CSV')),
                               tabPanel("Help / About",  htmlOutput('help')),
                             
                                if(admin){tabPanel("Custom Code",value="custom", aceEditor("ccrmd",wordWrap =T,mode="markdown",value="```{r}
2+2
plot(rnorm(23))
```"),   actionButton("cceval", "Run") ,htmlOutput("knitDoc"),downloadButton('report', 'Download PDF'), DT::dataTableOutput("ccfiles"),actionButton("ccload", "Load"),textInput("ccfilename",label="Filename"),actionButton("ccsave", "Save"),actionButton("ccrefresh", "Refresh"))}else{tabPanel(" ",  HTML("Nothing to see here"))}
                              
                  )))
  , verbatimTextOutput(outputId="debug") , plotOutput(outputId="debugplot"),HTML("<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-65298-21', 'auto');
  ga('send', 'pageview');

</script>")))
