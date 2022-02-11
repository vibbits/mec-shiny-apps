FROM rocker/r-ver:4.1.1

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/vibbits/nc-shiny-apps" \
      org.opencontainers.image.vendor="VIB Bioinformtics Core" \
      org.opencontainers.image.authors="Alexander Botzki <bits@vib.be>"

ENV S6_VERSION=v2.1.0.2
ENV SHINY_SERVER_VERSION=latest
ENV PANDOC_VERSION=default

ARG USER_ID=1000
ARG GROUP_ID=1000

RUN addgroup --gid $GROUP_ID shiny 
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID shiny 

COPY install_nc-shiny-apps.sh /rocker_scripts/
RUN ["chmod", "+x", "/rocker_scripts/install_nc-shiny-apps.sh"]

RUN /rocker_scripts/install_nc-shiny-apps.sh

USER shiny 

COPY --chown=shiny:shiny shiny-server/ /srv/shiny-server/

EXPOSE 3838

CMD ["/init"]

