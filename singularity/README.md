## Singularity RStudio Server

Image based on the Ubuntu 22.04 (R 4.3.2), see `./singularity/rstudio_server/recipe-4-3-2.def` for details

### First step

You need to clone this repo. Go to:

``` sh
$ cd singularity/
```

### Second step (Option 1): Download the image file:

https://zenodo.org/records/14251447

Note: There you will also find all the R packages used for the single cell analyses.

### Second step (Option 2): create the singularity image yourself:

``` sh
$ sudo singularity build rstudio-server.sif recipe-4-3-2.def
```

Note: The image must be created in `singularity/`, as it includes the file  `singularity/rstudio_auth.sh`

### Third step (config and run)

The following adjustments must be made in the `run.sh` file:

-   Set a port (variable `RSTUDIO_PORT`)
-   Set a path for the log files (variable `TMPDIR`)
-   The variable `RSTUDIO_PASSWORD` can be set to a password of your choice. Default is "password"
-   For parameter `auth-pam-helper` and `rsession-config-file` in run.sh: Set path to this folder.

Furthermore, set your R library path in `./singularity/rsession.conf`. I recommend to create a new folder. It can lead to conflicts when folders with libraries of another R version are used.

Now run the following :

``` sh
$ bash run.sh run rstudio-rstudio.sif
```

You can reach RStudio Server via your webbrowser (e.g. `localhost:8072/auth-sign-in`). Use the port that you have specified in the variable `RSTUDIO_PORT`. You can log in with your current user name and password you set in `RSTUDIO_PASSWORD`.

