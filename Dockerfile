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
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
COPY . /work

# Make every script executable (a fresh checkout / copy may drop the bit).
RUN find /work -type f \( -name '*.sh' -o -name '*.py' \) -exec chmod +x {} +

# The binaries are compiled by runme.sh when the container runs, since they are
# compiled with -march=native for the CPU that runs them.
CMD ["./runme.sh"]
