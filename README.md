# scTypeTC: single cell sub-typing of epithelial cells in thyroid cancer with glmnet-lasso model

The deadliest anaplastic thyroid cancer (ATC) often transforms from indolent differentiated thyroid cancer, most commonly, pappilary thyroid cancer (PTC). Large-scale single cell transcriptome data (scRNAseq) enabled the discovery of different epithelial and cancer cell subtypes occurred in different stages of thyrpod cancer progression.   scTypeTC is a computational tool using multinomial-lasso model to classify single cells with epithelial origins in thyroid cancer into 4 major subtypes, including normal thyroid follicular cells (TFCs), stress-responsive PTC cells (PTCs), inflammatory ATC cells (iATCs) and mesenchymal ATC cells (mATCs), using high throughput scRNAseq data of normal thyroids, PTC tumors and ATC tumors. scTypeTC takes our study results as ground truth to train the model that can be used to predict new scRNAseq datasets.  This vignette illustrate how to use the train model in scTypeTC to predict new single epithelial transcriptomes of thyroid tumors. 

## Step 1: installation
Installing scTypeTC from GitHub
```{r setup}
library(glmnet)
library(caret)
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
scTypeTC.test <- scTypeTC(test_data = test_data,cutoff_chance= 0.5,sam.name="test")
```

After this step, scTypeTC aumatically save the prediction results, the bar plot and dot plot in my working directory.  I can also extract the prediction result from scTypeTC.test object.


## Step 4: navigate prediction results

Now let's look at the prediction results.

```{r, eval=TRUE}
head(scTypeTC.test)
     
```
The first 4 columns in the prediction results matrix are the chance to be predicted as TFC, PTC, iATC or mATC. The last column is the final predicted subtype with the confident cutoff. Please note that cells predicted as undefined are in the results. The Rows are cells.
scTypeTC also report a matrix including the prediction powers of predictors.


