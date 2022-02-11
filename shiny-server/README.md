[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![shiny-apps](pictures/shiny.png) - Shiny-Apps
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

Please refer to the accompanying **[wiki](https://github.com/Nucleomics-VIB/shiny-apps/wiki)** for examples and workflows.

Read the header of each App.R to find out which packages you need to install on your server to have the code do its job.

## Shiny-apps
*[[back-to-top](#top)]*  

Those additional tools belong on a Shiny server and will execute R code in a interactive manner (please refder to https://shiny.rstudio.com/ for info about Shiny).


### **RNASeqFiltering.shinyapp** 
*[[Shiny-apps](#shiny-apps)]*

![RNASeqFiltering](pictures/RNASeqFiltering.png)

The **[RNASeqFiltering.shinyapp](RNASeqFiltering)** app loads a *StatisticalResults.xlsx* Excel file provided by the Nucleomics Core and filters it on one or more contrasts with chosen FDR and logFC limits. It then saves the results to XLSX and text files for further use. A live version was posted to https://nucleomics-core.shinyapps.io/RNASeqFiltering/. A sample ZIP file is present in the 'Data' subfolder for your convenience [(sample.zip)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/RNASeqFiltering/Data/sample.zip). A bundle can be downloaded [(here)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/RNASeqFiltering/Data/RNASeqFiltering-bundle.zip) to install the tool on your computer (you will need RStudio and a few R packages to run it)

### **fpkm2heatmap.shinyapp** 
*[[Shiny-apps](#shiny-apps)]*

![fpkm2heatmap](pictures/fpkm2heatmap.png)

The **[fpkm2heatmap.shinyapp](fpkm2heatmap)** app loads a FPKM Excel file provided by the Nucleomics Core and a list of EnsEMBL gene IDs (signature) and produced a heatmap plot taht can be tuned to your needs. A live version was posted to https://nucleomics-core.shinyapps.io/fpkm2heatmap/. A sample ZIP file is present in the 'Data' subfolder for your convenience [(sample.zip)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/fpkm2heatmap/Data/sample.zip). A bundle can be downloaded [(here)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/fpkm2heatmap/Data/fpkm2heatmap-bundle.zip) to install the tool on your computer (you will need RStudio and a few R packages to run it)

### **RBioanalyzer.shinyapp** 
*[[Shiny-apps](#shiny-apps)]*

The **[RBioanalyzer.shinyapp](RBioanalyzer)** app loads 2 to 3 Bioanalyzer exported csv files and creates an overlay plot. A live version was posted to https://nucleomics-core.shinyapps.io/RBioanalyzer/. A sample ZIP file is present in the 'Data' subfolder for your convenience [(sample.zip)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/RBioanalyzer/Data/sample.zip). A bundle can be downloaded [(here)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/RBioanalyzer/Data/RBioanalyzer-bundle.zip) to install the tool on your computer (you will need RStudio and a few R packages to run it)

### **RFilterRNASeq.shinyapp** 
*[[Shiny-apps](#shiny-apps)]*

![RFilterRNASeq](pictures/RFilterRNASeq.png)

The **[RFilterRNASeq.shinyapp](RFilterRNASeq)** app loads a StatisticalResults.xlsx file obtained from the Core, filters each contrast based on user input, and creates a Venn plot and a count table. The Venn plot supports up to 5 contrasts and is not created beyond that. A live version was posted to https://nucleomics-core.shinyapps.io/RFilterRNASeq/. A sample excel file with 2000 gene rows is present in the 'Data' subfolder for your convenience [(StatisticalResults.xlsx)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/RFilterRNASeq/Data/StatisticalResults.xlsx). A bundle can be downloaded [(here)](https://github.com/Nucleomics-VIB/Shiny-apps/raw/master/RFilterRNASeq/Data/RFilterRNASeq-bundle.zip) to install the tool on your computer (you will need RStudio and a few R packages to run it)


*[[back-to-top](#top)]*  

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
# nc-shiny-apps
