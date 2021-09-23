
3' Untranslated Regions of Tumor Suppressor Genes Evolved Specific Features to Favor Cancer Resistance
---------------------------------------------------------------------------------------

This resource provides the R code and processed data to reproduce key results described in Dan H, et al. 3' Untranslated Regions of Tumor Suppressor Genes Evolved Specific Features to Favor Cancer Resistance.

### Getting started
**1.** Clone Github repository. 
```
git clone https://github.com/Hd0909/Evolution-of-3UTR.git
```
**2.** Install required R packages.
```
Rscript install_required_packages.R
```
All the data required to generate the figures are in the ./data folder 
The bash script/Rscript  used to download or preprocess the raw data are also in the ./data folder


**3.** Run the following script to get the figures and supplementary figures
```
cd ./Rscript
#Supplementaryfigure1
Rscript Rscript_for_Supplementaryfigure1.R
## Figure 1
Rscript Rscript_for_figure1.R
# Figure 2,3 and supplementary Figure3-9
Rscript Rscript_for_figure23.R
# Figure 4,5 and supplementary Figure2,10-11
Rscript Rscript_for_figure45.R
# Figure 6 and supplementary Figure2,12-14
Rscript Rscript_for_figure6.R
# supplementary Figure15-16
Rscript Rscript_for_Supplementaryfigure15_16.R

```
### Contact
danhuang2018dana@gmail.com
