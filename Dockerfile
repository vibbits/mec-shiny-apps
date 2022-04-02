FROM rocker/r-ver:4.1.2

LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/vibbits/mec-shiny-apps" \
      org.opencontainers.image.vendor="VIB Bioinformtics Core" \
      org.opencontainers.image.authors="Alexander Botzki <bits@vib.be>"

ENV S6_VERSION=v2.1.0.2
ENV SHINY_SERVER_VERSION=latest
ENV PANDOC_VERSION=default

ARG USER_ID=1000
ARG GROUP_ID=1000

RUN addgroup --gid $GROUP_ID shiny 
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID shiny 

COPY install_mec-shiny-apps.sh /rocker_scripts/
RUN ["chmod", "+x", "/rocker_scripts/install_mec-shiny-apps.sh"]
RUN /rocker_scripts/install_mec-shiny-apps.sh

COPY register_fonts.sh /rocker_scripts/
RUN ["chmod", "ugo+x", "/rocker_scripts/register_fonts.sh"]

COPY register_fonts.r /usr/local/bin/
RUN ["chmod", "+x", "/usr/local/bin/register_fonts.r"]

USER shiny 

COPY --chown=shiny:shiny shiny-server/ /srv/shiny-server/
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN /rocker_scripts/register_fonts.sh

EXPOSE 3838

CMD ["/init"]
