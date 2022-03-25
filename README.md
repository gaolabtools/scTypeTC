# scTypeTC: Predict the subtype of thyroid epithelial cell and tumor cell into thyroid follicular cell, PTC, iATC and mATC

The deadliest anaplastic thyroid cancer (ATC) often transforms from indolent differentiated thyroid cancer.  We investigated this intra-tumor evolution process through integrative analysis of single cell transcriptomes and genetic alterations of thyroid cancer patients. We combined unsupervised clustering and ‘GLMnet’ machine learning to classify ATC tumor cells into two major subtypes: iATC and mATC.  The iATC cells over-expressed neutrophil-related inflammation genes and pathways. The mATC cells were characterized by gain of mesenchymal and anueploidy  phenotypes, and enrichment of ECM remodeling pathway.  scTypeTC is a computational tool using multinomial-lasso model to predict subtype of thyroid epithelial cell in single cell into thyroid follicular cell, PTC, iATC and mATC,following the best vote principles using sc-RNAseq data.  In this vignette, we will go through examples of how to use the scTypeTC R package to predict the sutype of single thyroid epithelial cell. 

## Step 1: installation
Installing scTypeTC from GitHub
```{r setup}
library(glmnet)
library(devtools)
install_github("gaolabtools/scTypeTC")
```

## Step 2: prepare normalized gene expression matrix input
The only one direct input that you need to prepare to run scTypeTC is the normalized gene expression matrix, with gene ids in rows and cell names in columns. The gene ids should be gene symbol. The cells should be thyroid epithelial cells. The matrix values can be easily generated from Seurat object. Below I provide an example of generating this normalized gene expression matrix from Seurat object.

### An example to generate input from Seurat object.
```{r, eval=FALSE}
  library(Seurat)
  test_data <- as.matrix(epi@assays$RNA@data)  #epi is a Seurat object of the thyroid epithelial cells
```

In this vignette, I take the example normalized gene expression matrix, test_data to demonstrate the workflow.

## Step 3: run scTypeTC
Now I have prepared the only one input, normalized gene expression matrix, I am ready to run scTypeTC. To filter out the prediction results with low confidence, I can tune parameters to keep only results with high chance of consistent predictions.  I put default cutoff=0.5. And cells below the cutoff are predicted as undefined. I can tune down this cutoff to keep more results in the prediction. I also give a sample name by setting sam.name="test". 

Now I runn the code:

```{r, message=FALSE}
library(scTypeTC)
scTypeTC.test <- scTypeTC(test_data = test_data,cutoff= 0.5,sam.name="test")
```

After this step, scTypeTC aumatically save the prediction results, the bar plot and dot plot in my working directory.  I can also extract the prediction result from scTypeTC.test object.


## Step 4: navigate prediction results

Now let's look at the prediction results.

```{r, eval=TRUE}
head(scTypeTC.test)
     
```
The first 4 columns in the prediction results matrix are the chance to be predicted as TFC, PTC, iATC or mATC. The last column is the final predicted subtype with the confident cutoff. Please note that cells predicted as undefined are in the results. The Rows are cells.
scTypeTC also generate a bar plot and a dot plot for the final predicted results, excluding the undefined cells.


