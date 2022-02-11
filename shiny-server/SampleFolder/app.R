# RFilterRNASeq.shinyapp
# A R/shiny tool to filter RNASeq data and create venn diagrams
#
# Stephane Plaisance, VIB Nucleomics Core
# visit our Git: https://github.com/Nucleomics-VIB
# version: 2017-06-14_v1.0
# © by using this tool, you accept the licence saved under ./www/licence.pdf

library("shiny")
library("openxlsx")
library("VennDiagram")
library("grDevices")

# you may uncomment the next line to allow large input files
options(shiny.maxRequestSize=1000*1024^2)
# the following test checks if we are running on shinnyapps.io to limit file size dynamically
# ref: https://stackoverflow.com/questions/31423144/how-to-know-if-the-app-is-running-at-local-or-on-server-r-shiny/31425801#31425801
#if ( Sys.getenv('SHINY_PORT') == "" ) { options(shiny.maxRequestSize=1000*1024^2) }

script.version="1.1"

# defaults variables and controllers
def.min.lfc <- 1
def.max.pv <- 0.05
chr.col <- NULL
all.contrasts <- NULL

# Define UI for application that draws a histogram
ui <- fluidPage(
  HTML('<style type="text/css">
    .row-fluid { width: 25%; }  
       .well { background-color: #99CCFF; }
       .shiny-html-output { font-size: 14px; line-height: 15px; }
       </style>'),
  # Application header
  headerPanel("Filter RNASeq data and create Venn plots"),

  # Application title
  titlePanel(
    windowTitle = "Filter RNASeq data and create Venn plots",
    tags$a(href="https://corefacilities.vib.be/nc", target="_blank",
           img(src='logo.png', align = "right", width="150", height="58.5", alt="VIB Nucleomics Core"))
  ),

  # Sidebar with input
  sidebarLayout(
    # show file import and molecule filters
    sidebarPanel(
      tags$h4(paste("code version: ", script.version, sep="")),
      downloadButton("downloadData", label = "Download test data"),
      tags$br(),
      tags$a(href="license.pdf", target="_blank", "usage licence"),
      tags$hr(),
      fileInput('file1', 'Choose XLSX File', accept='.xlsx'),
      tags$h4("modify cutoffs and click ", tags$em("Filter")),
      textInput('min.lfc', "minimal Fold-change (absolute value):", value = def.min.lfc),
      textInput('max.pv', "max corrected Pvalue (0..1):", value = def.max.pv),
      actionButton(inputId='goButton', "Filter", style='padding:4px; font-weight: bold; font-size:150%')
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput('plot1', width = "100%"),
      textOutput('full.data'),
      textOutput('min.lfc'),
      textOutput('max.pv'),
      textOutput('filt.data'),
      br(),
      tableOutput('filt.cnt')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$downloadData <- downloadHandler(
    filename <- function() { paste("StatisticalResults", "xlsx", sep=".") },
    content <- function(file) { file.copy("Data/StatisticalResults.xlsx", file) },
    contentType = "application/zip"
  )

  output$min.lfc <- renderText({
    paste("min log-FC (abs-val): ", input$min.lfc)
  })

  output$max.pv <- renderText({
    paste("max corrected-Pvalue: ", input$max.pv)
  })

  load.data <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)

    # load data from excel file
    dat <- read.xlsx(inFile$datapath, sheet=1)
    # count contasts
    chr.col <- which(colnames(dat)==as.vector("Chromosome"))
    # keep only filtering columns
    data <- dat[,c(2, 1, sort(c(seq(3, chr.col-1, 5),seq(5, chr.col-1, 5))))]
    # return data as 'load.data()'
    data
  })

  output$full.data <- reactive({
    if (is.null(load.data())) return(NULL)
    paste("gene rows in the Full data: ", nrow(load.data()))
  })

  filter.data <- eventReactive(input$goButton, {
    if (is.null(load.data())) return(NULL)
    min.lfc <- input$min.lfc
    max.fdr <- input$max.pv
    # store filtering results for each contrast
    num.contrasts <- (ncol(load.data()) - 2) / 2
    # indices for LR and FDR columns
    LR <- seq(1, 2 * num.contrasts, 2)
    FDR <- seq(2, 2 * num.contrasts, 2)
    # extract contrasts names
    names <- gsub(".:.logFC", "", colnames(load.data())[LR+2])
    # create list to store filtering results
    filtered.list <- vector("list", num.contrasts)
    # filter each contrast using cutoffs
    for (i in 1:num.contrasts){
      a <- LR[[i]]
      b <- FDR[[i]]
      # take one contrast and make small data.frame
      ct <- load.data()[,c(1, 2, a+2, b+2)]
      colnames(ct) <- c("Gene.Name", "Gene.ID", "LR", "FDR")

      filt <- ct[(abs(ct$LR)>=min.lfc & ct$FDR<=max.fdr),]
      filtered.list[[i]] <- as.vector(filt$Gene.ID)
    }
    # name contrasts in list
    names(filtered.list) <- names
    # return object
    filtered.list
  })

  output$filt.data <- reactive({
    if (is.null(filter.data())) return(NULL)
    paste("gene rows in the filtered data: ", length(unique(do.call(c, filter.data()))))
  })

  sum.data <- reactive({
    if (is.null(filter.data())) return(NULL)
    counts <- data.frame(do.call(cbind, lapply(filter.data(), length)),
                            row.names = "filtered gene counts")
    colnames(counts) <- names(filter.data())
    t(counts)
  })

  output$filt.cnt <- renderTable({sum.data()}, rownames=TRUE, colnames=TRUE, digits=2)

  output$plot1 <- renderPlot({
    if (is.null(filter.data())) return(NULL)
    if (length(filter.data()) > 5) return(NULL)
    # suppress log file creation
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    venn.plot <- venn.diagram(filter.data(),
                              filename = NULL,
                              names = names,
                              scaled = TRUE,
                              height = 800,
                              width = 800,
                              cex = 0.75,
                              cat.cex = 0.9,
                              margin=0.2,
                              cat.dist=rep(0.25, length(filter.data())))
    # plot
    grid.draw(venn.plot)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
