FROM rocker/hadleyverse

COPY r-packages.R /etc/r-packages.R
RUN R --slave < /etc/r-packages.R