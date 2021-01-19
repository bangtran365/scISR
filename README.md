# scISR: Single-cell Imputation using Subspace Regression
scISR performs imputation for single-cell sequencing data. scISR identifies the true dropout values in the scRNA-seq dataset using hyper-geomtric testing approach. Based on the result obtained from hyper-geometric testing, the original dataset is segregated into two including training data and imputable data. Next, training data is used for constructing a generalize linear regression model that is used for imputation on the imputable data.  
# How to install  
- The package can be installed from this repository.  
- Install devtools: `utils::install.packages('devtools')`  
- Install the package using: `devtools::install_github('bangtran365/scISR')`  
# Example   
- Load the sample dataset scISRExample: `data(scISRExample)`  
- Perform the : `imputed <- scISR(data = scISRExample$dropout)`  
- Plot the complete data: `plotscISR(scISRExample$raw, label = scISRExample$celltype)`  
- Plot the dropout data: `plotscISR(scISRExample$dropout, label = scISRExample$celltype)`  
- Plot the imputed data: `plotscISR(imputed, label = scISRExample$celltype)`  
