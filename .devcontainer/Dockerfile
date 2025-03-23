# Use official Python slim image for minimal size and better performance
FROM python:3.11-slim



# Install system dependencies required for building packages
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    procps \
    libssl-dev \
    zlib1g-dev \
    libffi-dev \
    libcurl4-openssl-dev \
    pkg-config && \
    rm -rf /var/lib/apt/lists/*

# Copy your requirements.txt
COPY requirements.txt cell_type_ontology.json /opt/



# Install Python packages
RUN pip install --upgrade pip && \
    pip install -r /opt/requirements.txt



# Default command (starts bash shell)
CMD ["bash"]