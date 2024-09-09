# Cell type gating models for scGate

Models can be loaded into scGate using the syntax:

``` r
models.db <- get_scGateDB()
model.Tcell <- models.db$human$generic$Tcell

gated <- scGate(obj, model=model.Tcell)
```

To use scGate as a multi-classifier, load one of the collections from the DB:

``` r
models.db <- get_scGateDB()
models.TME <- models.db$human$TME_HiRes

gated <- scGate(obj, model=models.TME)
table(gated$scGate_multi)
```

Note: we will tag versions of this DB using the format 'vN.M' (e.g. v1.4). This format is understood by scGate to download specific versions of the DB

## Model collections

For common use cases, there are ready-made model collections available. For example, there is a model collection to classify immune cells in the tumor micro-environment at high resolution:

``` r
scGate_models_DB <- get_scGateDB()
scGate_models <- scGate_models_DB$human$TME_HiRes
```

Model collections available:

-   **HiTME** <br> Default model used in the HiTME package. Includes immune cells at high resolution and stromal cells (e.g. fibroblasts).

-   **PBMC** <br> Includes immune cells at high resolution, optimized for possible background contamination from lysed blood cells.

-   **TME_broad** <br> Classification of stromal cells and immune cells at a broad level.

Changelog for model collections:
-   **PBMC**
  - For all cells, monocytes negative levels replaced with MoMacDC negative levels due to better specificity
  - Monocytes: Removed NK negative level because some Monocytes also express NKG7
  - Monocytes: kept only "myeloid" markers, as SPI1 is most specific for Monocytes