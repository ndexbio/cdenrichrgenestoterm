FROM continuumio/miniconda3

RUN pip install pandas
RUN pip install gseapy==0.10.4 --ignore-installed

RUN mkdir /tmp/cdenrichr
COPY dist/*.whl /tmp/cdenrichr/

RUN pip install /tmp/cdenrichr/cdenrichr*whl

RUN rm -rf /tmp/cdenrichr

ENTRYPOINT ["/opt/conda/bin/cdenrichrgenestoterm.py"]
CMD ["--help"]
