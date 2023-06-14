# Set base image
FROM tensorflow/tensorflow:2.1.0

LABEL program="PYTHIA"
LABEL description="A Deep Learning based Protein Blocks prediction tool."
LABEL version="1.2"
LABEL maintainer="gabriel.cretin@u-paris.fr"

# Set the working directory in the container
WORKDIR /PYTHIA

# Copy the dependencies file to the working directory
COPY requirements.txt .

# Install dependencies
RUN pip install --upgrade pip \
    && pip install -r requirements.txt

# Install the main script and models
COPY launch_PYTHIA.sh .
COPY models ./models

# Command to run on container start
ENTRYPOINT ["./launch_PYTHIA.sh"]
CMD ["--help"]
