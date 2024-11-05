docker build -t rmillikin/rnatools:latest . > docker_build.log 2>&1

# Push to dockerhub
docker push rmillikin/rnatools:latest