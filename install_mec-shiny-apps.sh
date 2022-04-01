#!/bin/bash
set -e

SHINY_SERVER_VERSION=${1:-${SHINY_SERVER_VERSION:-latest}}

# Run dependency scripts
. /rocker_scripts/install_s6init.sh
. /rocker_scripts/install_pandoc.sh

if [ "$SHINY_SERVER_VERSION" = "latest" ]; then
  SHINY_SERVER_VERSION=$(wget -qO- https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION)
fi

# Get apt packages
apt-get update
apt-get install -y --no-install-recommends \
    sudo \
    gdebi-core \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    make \
    zlib1g-dev
# <<<<<<< HEAD
#     zlib1g-dev \
#     make 
# =======
    

# >>>>>>> 249bc36881b147f8a3f4119661e56b541a61afa3

# Install Shiny server
wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-${SHINY_SERVER_VERSION}-amd64.deb" -O ss-latest.deb
gdebi -n ss-latest.deb
rm ss-latest.deb

# Get R packages
install2.r --error --skipinstalled -r "https://packagemanager.rstudio.com/cran/__linux__/focal/2022-01-28" shiny rmarkdown
# <<<<<<< HEAD
#install2.r --error --skipinstalled shiny rmarkdown
# =======
# >>>>>>> 249bc36881b147f8a3f4119661e56b541a61afa3

# Set up directories and permissions
#if [ -x "$(command -v rstudio-server)" ]; then
#  DEFAULT_USER=${DEFAULT_USER:-rstudio}
#  adduser ${DEFAULT_USER} shiny
#fi

cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/
chown shiny:shiny /var/lib/shiny-server
mkdir -p /var/log/shiny-server
chown shiny:shiny /var/log/shiny-server

# Doesn't work: Make standard linux truefonts available to shiny user
# chown shiny:shiny /usr/share/fonts/truetype

# Make truefonts available to shiny user as a copy
mkdir -p /var/fonts/truetype
cp -R /usr/share/fonts/truetype/* /var/fonts/truetype
chown shiny:shiny /var/fonts/truetype

# create init scripts
mkdir -p /etc/services.d/shiny-server
cat > /etc/services.d/shiny-server/run << 'EOF'
#!/usr/bin/with-contenv bash
## load /etc/environment vars first:
for line in $( cat /etc/environment ) ; do export $line > /dev/null; done
if [ "$APPLICATION_LOGS_TO_STDOUT" != "false" ]; then
    exec xtail /var/log/shiny-server/ &
fi
exec shiny-server 2>&1
EOF
chmod +x /etc/services.d/shiny-server/run

# Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## build ARGs
##NCPUS=${NCPUS:-1}

##apt-get update -qq && apt-get -y --no-install-recommends install \
##    libxml2-dev \
##    libcairo2-dev \
##    libgit2-dev \
##    default-libmysqlclient-dev \
##    libpq-dev \
##    libsasl2-dev \
##    libsqlite3-dev \
##    libssh2-1-dev \
##    libxtst6 \
##    libcurl4-openssl-dev \
##    unixodbc-dev && \
##  rm -rf /var/lib/apt/lists/*

#r packages
install2.r --error --skipinstalled -r "https://packagemanager.rstudio.com/cran/__linux__/focal/2020-06-04" \
    Rttf2pt1

install2.r --error --skipinstalled -r "https://packagemanager.rstudio.com/cran/__linux__/focal/2022-01-28" \
    remotes \
    shinyFeedback \
    shinyjs \
    shinyFiles \
    DT \
    colourpicker \
    here \
    vroom \
    forcats \
    dplyr \
    readr \
    tidyr \
    extrafont \
    extrafontdb \
    ggplot2 

# Make standard linux truefonts available to shiny user
# chown shiny:shiny /usr/local/lib/R/site-library/extrafontdb
# >>>>>>> 249bc36881b147f8a3f4119661e56b541a61afa3
## a bridge to far? -- brings in another 60 packages
# install2.r --error --skipinstalled -n $NCPUS tidymodels

 rm -rf /tmp/downloaded_packages
