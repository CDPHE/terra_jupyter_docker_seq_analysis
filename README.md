# terra-jupyter-python image
This repo contians a dockerfile for the docker image used for jupyter notebooks in the terra workspace.

This repo contains the terra-jupyter-python image that is compatible with notebook service in [Terra]("https://app.terra.bio/") called Leonardo. For example, use `us.gcr.io/broad-dsp-gcr-public/terra-jupyter-python:{version}` in terra. the dockerfile was origianlly cloned from https://github.com/DataBiosphere/terra-docker.git

## Image contents

The terra-jupyter-python image extends the [terra-jupyter-base](../terra-jupyter-base) image, supported in this repo.

To see the complete contents of this image please see the [Dockerfile](./Dockerfile).

I have appended iqtree (for maximum likelihood phylogenetic analysis) and the python module ete3 (for phylogenetic tree visualization) to the base python dockerfile.

See https://support.terra.bio/hc/en-us/articles/360037143432 for more info on creating a custom jupyter notebook environment for terra workspace.
