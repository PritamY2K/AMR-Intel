FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /app
COPY app_cloud.py /app/app_cloud.py
RUN mkdir -p /app/uploads /app/results /app/card_database

# Preload AMR databases during image build. This makes runtime faster and avoids first-run setup delays.
RUN cd /app/card_database && \
    wget -O card_data.tar.bz2 https://card.mcmaster.ca/latest/data && \
    tar -xjf card_data.tar.bz2 && \
    CARDJSON=$(find /app/card_database -name card.json | head -n 1) && \
    rgi load --card_json "$CARDJSON" --local && \
    amrfinder -u && \
    abricate --setupdb

EXPOSE 8501
ENV PATH=/opt/conda/bin:$PATH
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

CMD ["streamlit", "run", "app_cloud.py", "--server.address=0.0.0.0", "--server.port=8501"]
