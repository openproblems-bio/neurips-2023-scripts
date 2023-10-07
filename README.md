# neurips-2023-scripts
Scripts associated with the 2023 Open Problems competition

[![Run in Saturn Cloud](https://saturncloud.io/images/embed/run-in-saturn-cloud.svg)](https://app.community.saturnenterprise.io/dash/o/community/resources?templateId=312dc20b78914ed98f4b125fc70e00e8)

In order to run the notebook - you will probably need access to a GPU machine with 64GB of RAM. If you message Saturn Cloud (look for chat in the lower right window) and say "I'm competing in NeurIPS 2023", your account will be upgraded to allow access to more powerful instance types. Once you've been upgraded, you should be able to modify the instance tier to a machine with more RAM.

The Python image was built with the following environment.yaml

```
name: saturn
channels:
    - pytorch
    - fastai
    - rapidsai
    - nodefaults
    - conda-forge
dependencies:
    - blas=2.114=mkl
    - bokeh=2.4.3
    - cudatoolkit=11.3
    - dask-cuda=22.04
    - dask=2022.3.0
    - fastai=2.6.3
    - fsspec=2022.3.0
    - ipykernel=6.13.0
    - ipywidgets=7.7.0
    - matplotlib=3.5.2
    - numpy=1.21.6
    - pandas=1.4.2
    - pillow=9.1.1
    - pip
    - prefect=0.15.13
    - py-opencv=4.5.5
    - pyarrow=6.0.1
    - python-graphviz=0.20
    - python=3.9
    - pytorch=1.11.0
    - s3fs=2022.3.0
    - setuptools<60.0
    - tensorboard=2.9.0
    - torchaudio=0.11.0
    - torchvision=0.12.0
    - pynvml==11.4.1
    - pip:
          - dask-saturn>=0.4.1
          - prefect-saturn>=0.6.0
          - snowflake-connector-python==2.7.7
```

At startup, the template also installs the following packages with pip
```
pip install tables seaborn scanpy muon hdf5plugin plotly anndata kaggle scipy scikit-learn pandas matplotlib==3.7.2
```

There is also an separate R environment in the image that is used to run limma. R is installed via conda using the following environment.yaml

```
dependencies:
  - r-base=4.3.1=h29c4799_3
  - r-base64enc=0.1_3=r43h57805ef_1006
  - r-cli=3.6.1=r43ha503ecb_1
  - r-cpp11=0.4.6=r43hc72bb7e_0
  - r-crayon=1.5.2=r43hc72bb7e_2
  - r-digest=0.6.33=r43ha503ecb_0
  - r-ellipsis=0.3.2=r43h57805ef_2
  - r-evaluate=0.21=r43hc72bb7e_1
  - r-fansi=1.0.4=r43h57805ef_1
  - r-fastmap=1.1.1=r43ha503ecb_1
  - r-glue=1.6.2=r43h57805ef_2
  - r-htmltools=0.5.6=r43ha503ecb_0
  - r-irdisplay=1.1=r43hd8ed1ab_2
  - r-irkernel=1.3.2=r43h785f33e_1
  - r-jsonlite=1.8.7=r43h57805ef_0
  - r-lifecycle=1.0.3=r43hc72bb7e_2
  - r-pbdzmq=0.3_9=r43ha81a24b_1
  - r-pillar=1.9.0=r43hc72bb7e_1
  - r-ragg=1.2.5=r43h85cdef0_2
  - r-repr=1.1.6=r43h785f33e_1
  - r-rlang=1.1.1=r43ha503ecb_1
  - r-systemfonts=1.0.4=r43haf97adc_2
  - r-textshaping=0.3.6=r43h24cd192_6
  - r-utf8=1.2.3=r43h57805ef_1
  - r-uuid=1.1_1=r43h57805ef_0
  - r-vctrs=0.6.3=r43ha503ecb_0
```

After that, an R script is used to install the final dependencies

```
install.packages("remotes", repos = "http://cran.us.r-project.org")
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("edgeR")
remotes::install_cran(c("anndata"), repos = "https://cran.rstudio.com")
install.packages("optparse", repos = "http://cran.us.r-project.org")
remotes::install_version("reticulate", version = "1.22", repos = "http://cran.us.r-project.org")
```
