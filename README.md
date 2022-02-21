# NeoEnrichment
R package to do enrichment analysis for neoantigens
Install by `devtools::install_github("wt12318/NeoEnrichment",ref="dev")`

## Usage

For calculating ESccf, we need supply a dataframe with at least 3 columns:

- sample, sample name 
- neo, indicating whether a mutation (a row of the dataframe) is neoantigentic mutation; value can be "yes" or "no"
- ccf, indicating the CCF value of mutations, range from 0 to 1.

Then we can use the `cal_nes_new_test` function to calculate ESccf for one sample:

```{r}
a <- NeoEnrichment::cal_nes_new_test(dt = data, sample_counts = 1000, need_p = FALSE)
```

There are three parameters of the function `cal_nes_new_test`:

- dt, the mutation dataframe mentioned above
- need_p, whether need calculated p values
- sample_counts, the number of random sampling when calculate p values.

For calculating ESrna, we also need supply a dataframe with at least 3 columns:
- sample, sample name 
- neo, indicating whether a mutation (a row of the dataframe) is neoantigentic mutation; value can be "neo" or "not_neo"
- exp, indicating the expression of gene which the mutation located (often in TPM unit)

Then we can use the `cales_t` function to calculate ESexp for samples (can be used for multiple samples):

```{r}
a <- NeoEnrichment::cales_t(data = dt,barcode = x,type = "II",
                            calp = FALSE,sample_counts = 1000,
                            cal_type = "exp")
```

There are six parameters of the function `cales_t`:

- data, the mutation dataframe mentioned above
- barcode, the barcode of the sample needed to run
- type, "I" or "II", "I" means put more weight on neoantigentic mutations, while "II" means put equal weights, we used "II" in our paper
- calp, whether need calculated p values,
- sample_counts, the number of random sampling when calculate p values
- cal_type, the ES type we need calculate, the "CCF" was discarded.




