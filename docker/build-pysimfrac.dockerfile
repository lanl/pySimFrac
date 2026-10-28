##=====================================================##
## pySimFrac + LaGriT Dockerfile (multi-stage, multi-arch)
##=====================================================##

# ---------- Stage 1: build LaGriT ----------
FROM ubuntu:22.04 AS lagrit-builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        git \
        m4 \
        bison \
        zlib1g-dev \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
RUN git clone --depth 1 https://github.com/lanl/LaGriT.git && \
    # Strip x86-only flags so the build works on arm64 too.
    # -m64 is an x86_64-vs-x86_32 selector that gcc/gfortran reject on arm64.
    # -Dlinx64 / -Dmacx64 are LaGriT internal arch defines that we keep on x86 only.
    find LaGriT -type f \( -name 'CMakeLists.txt' -o -name '*.cmake' \) \
        -exec sed -i 's/-m64//g; s/-Dlinx64//g; s/-Dmacx64//g' {} + && \
    # On x86_64, put -Dlinx64 back (LaGriT source uses it for arch-specific code paths)
    if [ "$(uname -m)" = "x86_64" ]; then \
        find LaGriT -type f \( -name 'CMakeLists.txt' -o -name '*.cmake' \) \
            -exec sed -i 's/-fcray-pointer/-Dlinx64 -fcray-pointer/' {} + ; \
    fi && \
    cmake -S LaGriT -B LaGriT/build && \
    cmake --build LaGriT/build -j"$(nproc)"


# ---------- Stage 2: runtime image ----------
FROM ubuntu:22.04 AS runtime

LABEL maintainer="Daniel Livingston <livingston@lanl.gov>"

ENV DEBIAN_FRONTEND=noninteractive
ENV APP_PATH=/lagrit
ENV PATH=$APP_PATH/bin:$PATH

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-tk \
        libgfortran5 \
        libgomp1 \
        libstdc++6 \
        zlib1g \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p $APP_PATH/bin
COPY --from=lagrit-builder /build/LaGriT/build/lagrit $APP_PATH/bin/lagrit

WORKDIR /app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . ./
ENV PYTHONPATH=/app/src

CMD ["bash"]