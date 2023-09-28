FROM        ubuntu:20.04
ENV         LC_ALL=C TZ=Asia/Shanghai PS=20201110 PT=9026
RUN         ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
            apt update && apt upgrade -y && apt install openssh-server language-pack-zh-hans language-pack-zh-hans-base aria2 -y && apt autoremove -y
RUN         sed -i "s/#PermitEmptyPasswords no/PermitEmptyPasswords no/g" /etc/ssh/sshd_config && \
            sed -i "s/#PermitRootLogin prohibit-password/PermitRootLogin yes/g" /etc/ssh/sshd_config && \
            sed -i "s/#ListenAddress 0.0.0.0/ListenAddress 0.0.0.0/g" /etc/ssh/sshd_config && \
            sed -i "s/#LoginGraceTime 2m/LoginGraceTime 2m/g" /etc/ssh/sshd_config && \
            echo root:${PS} | chpasswd
COPY        --chown=root:root ./Mambaforge-pypy3-Linux-x86_64.sh /root/Mambaforge-pypy3-Linux-x86_64.sh
RUN         bash  ~/Mambaforge-pypy3-Linux-x86_64.sh -b && \
            rm -f ~/Mambaforge-pypy3-Linux-x86_64.sh    && \
            ~/mambaforge-pypy3/bin/mamba init
COPY        --chown=root:root ./condarc /root/.condarc
ENTRYPOINT  ["/bin/bash", "-c", "echo root:${PS} | chpasswd && sed -i \"s/#Port 22/Port ${PT}/g\" /etc/ssh/sshd_config && service ssh start -D"]

#docker build     -t    doujiangbaozi/sliverworkspace:1.06 .
#docker tag       bd51  doujiangbaozi/sliverworkspace:latest
#docker push            doujiangbaozi/sliverworkspace:1.05
#docker push            doujiangbaozi/sliverworkspace:latest

#EXPOSE      $port
#sed -i "s/#Port 22/Port ${port}/g" /etc/ssh/sshd_config && \

#~/mambaforge-pypy3/bin/mamba config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free         && \
#~/mambaforge-pypy3/bin/mamba config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge && \
#~/mambaforge-pypy3/bin/mamba config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda    && \
#~/mambaforge-pypy3/bin/mamba config --set show_channel_urls yes    && \

#aria2c https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-Linux-x86_64.sh -d .