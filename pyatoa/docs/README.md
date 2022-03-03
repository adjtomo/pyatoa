The Pyatoa docs are defined by a combination of reStructuredText (.rst) files
and Jupyter/IPython notebooks (.ipynb). These docs pages are hosted on 
ReadTheDocs (RTD). 

To reduce complexity server-side on RTD, the docs are built locally and then
uploaded to GitHub which triggers RTD to build them. Juypter notebooks are NOT
used on RTD, all files are expected as .rst files.

To build the .rst files, you must run the following command, which will execute
all code within Jupyter notebooks (\*.ipynb), generate intermediate .rst files. 
These files and the corresponding \*\_files/ directories, which contain .png 
images necessary for the docs, must then be uploaded to the public repository.

```bash
python convert.py
```

If you would like to bulid the documentation locally, you can run the following:

```bash
make html
```

which will build the docs into the \_build/ directory. You can then open
\_build/html/index.html to look at the docs.

These docs rely on sphinx-autoapi, which is called during the Make command. 
This autogenerates API without needing to install the package, keeping things
more lightweight. sphinx-autoapi will automatically include a link to the 
API reference in the index toctree.
