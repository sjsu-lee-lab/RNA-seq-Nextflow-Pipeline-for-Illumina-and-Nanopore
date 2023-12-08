#!/bin/bash

cd $1

# Nextflow requires  Java 11 (or later, up to 21) to be installed

# Install SDK to get Java
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"

# Install Java
sdk install java 17.0.6-amzn

# Install Nextflow
wget -qO- https://get.nextflow.io | bash