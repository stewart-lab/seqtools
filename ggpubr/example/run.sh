# wait 3 seconds and then open browser as a background process.
# this keeps the terminal free to run the docker command and the browser open.
# I wrote it this way so Docker is the foreground process and you can exit it with ctrl+c.

# Use the following credentials to log in:
#   username: rstudio
#   password: password
(sleep 3 && open http://localhost:8787) &

# Run the docker container
docker run --rm -v "$(pwd)/data":/mnt/data -p 8787:8787 -e PASSWORD=password rmillikin/ggpubr:latest