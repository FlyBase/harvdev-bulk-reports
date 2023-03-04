FROM flybase/harvdev-docker:latest

WORKDIR /src

RUN mkdir /src/input/
RUN mkdir /src/logs/
RUN mkdir /src/output/
RUN mkdir /src/temp/

ADD src/**              harvdev-bulk-reports/
ADD requirements.txt    harvdev-bulk-reports/requirements.txt
ADD perl_modules/**     /usr/local/lib/perl5/site_perl/

RUN pip3 install -r harvdev-bulk-reports/requirements.txt

ENTRYPOINT [ "/bin/bash" ]
