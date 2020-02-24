
# install.packages("stevedore")
library("stevedore")
docker <- stevedore::docker_client()
docker$container$run("nextflow/nextflow:latest", c("echo", "hello world"))  # , detach = TRUE)

docker$container$list()
# docker$container$

docker$image$list()
