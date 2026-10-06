FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        make \
        bash \
        wget \
        curl \
        gzip \
        gawk \
        tar \
        unzip \
        ca-certificates \
        default-jre-headless \
        gnuplot-nox \
        texlive-latex-base \
        texlive-pictures \
        texlive-latex-recommended \
        python3 \
        python3-pip \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir --break-system-packages \
        pandas matplotlib seaborn

WORKDIR /work
COPY . /work

# Make every script executable (a fresh checkout / copy may drop the bit).
RUN find /work -type f \( -name '*.sh' -o -name '*.py' \) -exec chmod +x {} +

RUN make -C /work/experiments/src

CMD ["./runme.sh"]
