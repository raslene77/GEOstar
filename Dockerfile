FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    cmake g++ gcc git make ninja-build pkg-config zip unzip curl \
    libglfw3-dev libglm-dev libgl1-mesa-dev

RUN git clone https://github.com/microsoft/vcpkg.git /opt/vcpkg && \
    /opt/vcpkg/bootstrap-vcpkg.sh && \
    /opt/vcpkg/vcpkg install glad

WORKDIR /app
COPY . .

RUN cmake -B build -S . \
    -DCMAKE_TOOLCHAIN_FILE=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake && \
    cmake --build build -j

CMD ["/bin/bash"]