# Pyatoa Documentation

The Pyatoa documentation is built with Sphinx and ReadTheDocs. The official 
documentation can be found at https://pyatoa.readthedocs.io. 

Docs building is automatically triggered when updates are pushed to Pyatoa. 
Each branch of the code may have a different documentation. The 'latest' 
version of the docs points to the 'devel' branch of the code.

## Building Docs Locally

In order to build the Docs locally, you will first need to create a separate 
Conda environment with a few packages, you can do this by running:

``` bash
conda env create --file environment.yaml
conda activate pyatoa-docs
```

You can then run the make command to generate the .html files. You can find your 
local docs in the *_build/html* directory

```bash
make html
```

## Docs from Jupyter Notebooks

Some docs pages are generated from Jupyter notebooks. Typically, these are docs 
which execute actual code blocks and display their output 
(e.g., tutorial walkthroughs). 

The corresponding .rst page are then converted from the Jupyter notebooks using
a Python script. The *notebooks/* directory contains docs pages that are 
generated in this way. To generate a .rst file from a Jupyter notebook

```bash
cd notebooks/
python convert.py <FILENAME>  # replace filename with any of the available .ipynb files
```

This automatically deletes the current .rst file and replaces it with the newly generated
file, converted from the notebook. Additional functions in the convert script
will edit some of the rst commands to make the resulting docs page look a bit 
cleaner.

> __NOTE:__ The notebook is **not** executed with the convert script, it only 
> takes the current state of the notebook and converts it to RST. If you want 
> to edit cell inputs or outputs you will need to install Pyatoa to the Conda 
> environment and execute the notebook on your own.
