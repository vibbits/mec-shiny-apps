#!/usr/bin/env r

# specific for MEX linux server, install extrafontdb in accessible folder
if (!dir.exists('/home/shiny/.fonts')) {
    dir.create('/home/shiny/.fonts')
    file.copy("/home/shiny/www/DejaVuSans.ttf", "/home/shiny/.fonts")
    system('fc-cache -f /home/shiny/.fonts')

    .libPaths(c('/srv/shiny-server/TraVisPies/r-lib', .libPaths()))
    install.packages('/srv/shiny-server/TraVisPies/r-lib/extrafontdb_1.0.tar.gz',type = 'source',repos = NULL)
}

library(Rttf2pt1)
library(extrafont)    #for using other fonts
font_import(prompt=F)
