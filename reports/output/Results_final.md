---
title: "Recombination Landscapes in Angiosperms - Final results"
date: "2022-06-21"
csl: elsevier-harvard2.csl
output:
  # bookdown::pdf_document2: null
  bookdown::html_document2:
    keep_md: true
    toc: yes
  # pdf_document: default
  # bookdown::pdf_document2:
  #   toc: yes
  # bookdown::word_document2: default
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding,
  output_dir = "output", output_format = "all") })
bibliography: bibliography.bib
---

<style>
body {
text-align: justify}
</style>




# Dataset








Save metadatas about the dataset.



#### Figure S1

Save Marey maps in a Supplementary pdf...









#### Figure S2

Save recombination maps in a Supplementary pdf...













## Dataset quality

Assessment of recombination landscapes quality.



```r
# ============================================================================
# Quantitative assessment - Correlation between the chromosome wide recombination
# rate & the mean estimated recombination rate Show the linear relationship
# between chromosome wide rate and the mean estimated recombination rate as a
# cross-validation procedure
# ============================================================================
plot(chromosomes$chrwide.rate, chromosomes$mean.recrate, xlab = "Chromosome-wide recombination rate", 
    ylab = "Mean estimated recombination rate", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2)
abline(a = 0, b = 1, col = "Red")
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-9-1.png" alt="(ref:landscape-quality)"  />
<p class="caption">(\#fig:unnamed-chunk-9)(ref:landscape-quality)</p>
</div>

(ref:landscape-quality) Correlation between the chromosome wide recombination rate and the mean estimated recombination rate. Identity line drawn in red.










```r
cor.test(chromosomes$chrwide.rate, chromosomes$mean.recrate, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  chromosomes$chrwide.rate and chromosomes$mean.recrate
## S = 117170, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9976094
```

Correlation between recombination maps at 100kb and 1Mb.




```r
# Correlation
cor.test(estimates1Mb$rec.rate, estimates1Mb$average100kb, method = "spearman")
```

```
## Warning in cor.test.default(estimates1Mb$rec.rate, estimates1Mb$average100kb, :
## Cannot compute exact p-value with ties
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  estimates1Mb$rec.rate and estimates1Mb$average100kb
## S = 9.5475e+10, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9975762
```

```r
plot(estimates1Mb$rec.rate, estimates1Mb$average100kb)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-13-1.png)<!-- -->




```r
sp_corr$species = gsub("_", " ", sp_corr$species)
colnames(sp_corr) = c("Species", "Number of 1Mb windows", "Rho", "Statistic", "p-value")
kable(sp_corr, digits = 3, align = "c")
```

<table>
 <thead>
  <tr>
   <th style="text-align:center;"> Species </th>
   <th style="text-align:center;"> Number of 1Mb windows </th>
   <th style="text-align:center;"> Rho </th>
   <th style="text-align:center;"> Statistic </th>
   <th style="text-align:center;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> 1566 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 38660.408 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> 121 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 80.001 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> 805 </td>
   <td style="text-align:center;"> 0.998 </td>
   <td style="text-align:center;"> 204003.923 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> 2018 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 315.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Boechera stricta </td>
   <td style="text-align:center;"> 189 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 634.535 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> 138 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 10.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> 382 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 2404.004 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> 179 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 40.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> 447 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 6661.284 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> 994 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 3069.505 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> 111 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 250.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsicum annuum </td>
   <td style="text-align:center;"> 2612 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 1116249.438 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cenchrus americanus </td>
   <td style="text-align:center;"> 1535 </td>
   <td style="text-align:center;"> 0.932 </td>
   <td style="text-align:center;"> 38673978.897 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrullus lanatus </td>
   <td style="text-align:center;"> 309 </td>
   <td style="text-align:center;"> 0.998 </td>
   <td style="text-align:center;"> 8987.501 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> 251 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 263.030 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Coffea canephora </td>
   <td style="text-align:center;"> 587 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 43672.649 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> 339 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 1637.193 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> 193 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 356.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> 201 </td>
   <td style="text-align:center;"> 0.985 </td>
   <td style="text-align:center;"> 14935.074 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> 182 </td>
   <td style="text-align:center;"> 0.972 </td>
   <td style="text-align:center;"> 20934.056 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> 343 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 1372.518 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Draba nivalis </td>
   <td style="text-align:center;"> 242 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 146.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> 519 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 4390.011 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> 614 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 2980.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> 922 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 66360.603 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> 1941 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 276988.789 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> 756 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 6256.510 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> 2240 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 935972.729 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> 4581 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 4321294.375 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Juglans regia </td>
   <td style="text-align:center;"> 721 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 13974.110 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> 442 </td>
   <td style="text-align:center;"> 0.998 </td>
   <td style="text-align:center;"> 28967.238 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> 564 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 24885.071 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> 656 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 1636.004 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Mangifera indica </td>
   <td style="text-align:center;"> 363 </td>
   <td style="text-align:center;"> 0.998 </td>
   <td style="text-align:center;"> 12414.386 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> 391 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 5484.362 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Momordica charantia </td>
   <td style="text-align:center;"> 297 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 3431.625 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Nelumbo nucifera </td>
   <td style="text-align:center;"> 812 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 4921.502 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> 373 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 761.005 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> 377 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 1280.500 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> 488 </td>
   <td style="text-align:center;"> 0.991 </td>
   <td style="text-align:center;"> 173277.055 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> 520 </td>
   <td style="text-align:center;"> 0.998 </td>
   <td style="text-align:center;"> 39531.933 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> 180 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 962.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> 227 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 96.500 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Quercus sp </td>
   <td style="text-align:center;"> 720 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 7249.008 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Raphanus sativus </td>
   <td style="text-align:center;"> 232 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 553.352 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> 259 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 973.511 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> 405 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 154.004 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> 806 </td>
   <td style="text-align:center;"> 0.995 </td>
   <td style="text-align:center;"> 374629.627 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> 728 </td>
   <td style="text-align:center;"> 0.999 </td>
   <td style="text-align:center;"> 82532.028 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> 663 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 7314.467 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> 335 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 272.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> 13346 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 39822683.941 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum dicoccoides </td>
   <td style="text-align:center;"> 5821 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 1611892.897 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum urartu </td>
   <td style="text-align:center;"> 4611 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 441062.644 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> 478 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 191.000 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> 467 </td>
   <td style="text-align:center;"> 0.988 </td>
   <td style="text-align:center;"> 182751.865 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> 2107 </td>
   <td style="text-align:center;"> 1.000 </td>
   <td style="text-align:center;"> 20666.063 </td>
   <td style="text-align:center;"> 0 </td>
  </tr>
</tbody>
</table>

```r
# Save the table of correlation per species/dataset
write.xlsx(x = sp_corr, file = paste(wd, "tables/article_one_supp/tables_supplementary.xls", 
    sep = ""), sheetName = "S4_MareyMaps_CorrelationScales", row.names = FALSE, append = TRUE)
```


### Correction of the linkage map length uncorrected/corrected (map coverage)

Two methods to correct linkage map length
1. Chakravarti et al. (1991) 
2. Hall and Willis (2005)


<!-- ```{r echo = TRUE} -->
<!-- # Two methods to correct linkage map length -->
<!-- # 1. Chakravarti et al. (1991) -->
<!-- cor.test(chromosome.stats$linkage.map.length, chromosome.stats$linkage.map.length.correctedCh, method = "spearman") -->
<!-- plot(chromosome.stats$linkage.map.length, chromosome.stats$linkage.map.length.correctedCh,  -->
<!--      xlab = "Linkage map length (uncorrected)", -->
<!--      ylab = "Linkage map length (corrected)") -->

<!-- #  2. Hall and Willis (2005) -->
<!-- cor.test(chromosome.stats$linkage.map.length, chromosome.stats$linkage.map.length.correctedHW, method = "spearman") -->
<!-- plot(chromosome.stats$linkage.map.length, chromosome.stats$linkage.map.length.correctedHW,  -->
<!--      xlab = "Linkage map length (uncorrected)", -->
<!--      ylab = "Linkage map length (corrected)") -->
<!-- ``` -->


Corrected linkkage map length don't change significantly the raw linkage map length.


```r
corr_linkagemap = chromosome.stats$linkage.map.length - chromosome.stats$linkage.map.length.correctedHW
summary(corr_linkagemap)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -5.6216 -1.7449 -0.8285 -1.1891 -0.3466  0.0000
```

```r
n.boot = 1000
boot = numeric(n.boot)
for (i in 1:n.boot) {
    boot[i] = mean(sample(corr_linkagemap, replace = TRUE), na.rm = TRUE)
}
mean(boot)
```

```
## [1] -1.184658
```

```r
quantile(boot, c(0.025, 0.975))
```

```
##      2.5%     97.5% 
## -1.271591 -1.096404
```


### Genome coverage



```r
# The difference between genomic map length (Mb) in the Marey map and the total
# chromosome size from the fasta file.
summary(chromosome.stats$phys.map.length * 10^6 - chromosome.stats$chromosome.size.bp)
```

```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -1.490e-08  0.000e+00  0.000e+00  3.081e-11  0.000e+00  7.451e-09
```



### Marker density

Number of markers


```r
summary(chromosome.stats$nb.markers)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    31.0   117.0   224.0   956.6   581.0 49483.0
```

```r
hist(chromosome.stats$nb.markers, breaks = 300, main = "", xlab = "Number of markers")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
hist(chromosome.stats$nb.markers, breaks = 30000, xlim = c(0, 1000), main = "", xlab = "Number of markers")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-18-2.png)<!-- -->


Marker density (nb of markers/total map length); i.e. the averaged number of markers per cM/Mb, more is better.


```r
summary(chromosome.stats$density.markers.cM)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.3571   1.0172   2.2403   8.3394   4.6830 359.1491
```

```r
hist(chromosome.stats$density.markers.cM, breaks = 80, main = "", xlab = "Marker density (cM)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-19-1.png)<!-- -->



```r
summary(chromosome.stats$density.markers.bp * 10^6)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.1563   2.4660   4.8309  20.3394  14.2718 696.2378
```

```r
hist(chromosome.stats$density.markers.bp * 10^6, breaks = 80, main = "", xlab = "Marker density (Mb)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
hist(chromosome.stats$density.markers.bp * 10^6, breaks = 800, main = "", xlab = "Marker density (Mb)", 
    xlim = c(0, 100))
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-20-2.png)<!-- -->

Mean interval between markers (cM/Mb)

TODO - median interval (less sensitive to extreme values, i.e. large gaps)


```r
summary(chromosome.stats$marker.interval.cM)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.002785 0.212289 0.447341 0.659125 0.995402 2.810800
```

```r
hist(chromosome.stats$marker.interval.cM, breaks = 80, main = "", xlab = "Mean interval between adjacent markers (cM)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-21-1.png)<!-- -->



```r
summary(chromosome.stats$marker.interval.bp/10^6)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.001436 0.069858 0.207525 0.341234 0.403838 6.522944
```

```r
hist(chromosome.stats$marker.interval.bp/10^6, breaks = 80, main = "", xlab = "Mean interval between adjacent markers (Mb)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
hist(chromosome.stats$marker.interval.bp/10^6, breaks = 80, main = "", xlab = "Mean interval between adjacent markers (Mb)", 
    xlim = c(0, 2))
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-22-2.png)<!-- -->



### How many maps discarded

Number of dataset discarded.


Number of chromosomes discarded among retained dataset.


```r
maps = read.table(file = paste0(wd, "/data-cleaned/marey_maps/AllMaps.txt"), header = TRUE, 
    sep = "\t")

# All chromosomes - import Marey maps
set_list = paste(wd, "/data-cleaned/marey_maps/", unique(maps$set), ".txt", sep = "")
allmaps = do.call("rbind", lapply(set_list, FUN = function(file) {
    read.table(file, header = TRUE, sep = "\t")
}))
# Chromosomes in dataset
length(unique(paste(allmaps$set, allmaps$map, sep = "_")))
```

```
## [1] 682
```

```r
# Chromosomes retained
length(unique(paste(maps$set[which(maps$vld)], maps$map[which(maps$vld)], sep = "_")))
```

```
## [1] 665
```




## Phylogeny

Phylogenetic tree of sampled species.


```
## 
## Phylogenetic tree with 57 tips and 56 internal nodes.
## 
## Tip labels:
##   Manihot esculenta, Cucumis sativus, Cucumis melo, Citrullus lanatus, Cucurbita pepo, Cucurbita maxima, ...
## Node labels:
##   Root, Pentapetalae, mrcaott2ott8384, rosids, mrcaott2ott371, mrcaott371ott579, ...
## 
## Rooted; includes branch lengths.
```

#### Figure S3

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-25-1.png" alt="(ref:phylogenetic-tree-annotated)"  />
<p class="caption">(\#fig:unnamed-chunk-25)(ref:phylogenetic-tree-annotated)</p>
</div>

(ref:phylogenetic-tree-annotated) Phylogenetic tree representing sampled species and annotated with mean recombination rates and mean chromosome sizes.






## Life history traits

### Mating system


```r
lht = openxlsx::read.xlsx(paste(wd, "data-cleaned/ListPlant_Thomas.xlsx", sep = ""), 
    sheet = 1)
lht$sp = paste(lht$Genus, " ", lht$Species, sep = "")
chr_lht = chromosome.stats
chr_lht$mating = NA
for (i in 1:nrow(chr_lht)) {
    chr_lht$mating[i] = lht$Mating_system[which(lht$sp == chr_lht$species[i])]
}
lht$recrate = NA
lht$linkagemaplength = NA
for (i in 1:nrow(lht)) {
    lht$recrate[i] = mean(chr_lht$mean.recrate[which(chr_lht$species == lht$sp[i])], 
        na.rm = TRUE)
    lht$linkagemaplength[i] = mean(chr_lht$linkage.map.length.correctedHW[which(chr_lht$species == 
        lht$sp[i])], na.rm = TRUE)
}
# Sample sizes
table(lht$Mating_system)
```

```
## 
##       Mixed Outcrossing     Selfing 
##          12          21          20
```

```r
table(chr_lht$mating)
```

```
## 
##       Mixed Outcrossing     Selfing 
##         163         200         211
```

```r
boxplot(recrate ~ Mating_system, data = lht, xlab = "Mean recombination rate", ylab = "Linkage map length (cM)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```r
kruskal.test(recrate ~ Mating_system, data = lht)
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  recrate by Mating_system
## Kruskal-Wallis chi-squared = 1.7842, df = 2, p-value = 0.4098
```

```r
boxplot(linkagemaplength ~ Mating_system, data = lht, xlab = "Mating system", ylab = "Linkage map length (cM)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-27-2.png)<!-- -->

```r
kruskal.test(linkagemaplength ~ Mating_system, data = lht)
```

```
## 
## 	Kruskal-Wallis rank sum test
## 
## data:  linkagemaplength by Mating_system
## Kruskal-Wallis chi-squared = 0.90034, df = 2, p-value = 0.6375
```

```r
boxplot(mean.recrate ~ mating, data = chr_lht, xlab = "Mean recombination rate", 
    ylab = "Linkage map length (cM)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-27-3.png)<!-- -->

```r
boxplot(linkage.map.length ~ mating, data = chr_lht, xlab = "Mating system", ylab = "Linkage map length (cM)")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-27-4.png)<!-- -->


## Position of the centromeric index on the genome


False Positive Rate. How many centromeric indexes have been assigned to the wrong side of the chromosome? Compare the recombination rate at the centromeric index with the putative position on the opposite direction, thus the FPR is the percentage of CI that have a higher recombination rate than the opposite position.






```r
sum(!is.na(df$true_CI_rate))
```

```
## [1] 424
```

```r
# Number of centromeric indexes wrongly oriented (i.e. lower recombination rate
# on the opposite position)
sum(df$true_CI_rate > df$false_CI_rate, na.rm = TRUE)
```

```
## [1] 70
```

```r
# Proportion -> the FPR
sum(df$true_CI_rate > df$false_CI_rate, na.rm = TRUE)/sum(!is.na(df$true_CI_rate))
```

```
## [1] 0.1650943
```

```r
# Identify these chromosomes to put them asides
idx_wrongCIorientation = which(df$true_CI_rate > df$false_CI_rate)
paste(df$set[which(df$true_CI_rate > df$false_CI_rate)], df$chromosome[which(df$true_CI_rate > 
    df$false_CI_rate)])
```

```
##  [1] "Boechera_stricta_Lee2020 2"                
##  [2] "Boechera_stricta_Lee2020 5"                
##  [3] "Brassica_napus_Yang2017 A07"               
##  [4] "Brassica_rapa_CodyMarkelz2017 A01"         
##  [5] "Capsella_rubella_Slotte2013 1"             
##  [6] "Capsella_rubella_Slotte2013 4"             
##  [7] "Capsicum_annuum_Han2016 3"                 
##  [8] "Capsicum_annuum_Han2016 12"                
##  [9] "Cenchrus_americanus_Pucher2017 2"          
## [10] "Cenchrus_americanus_Pucher2017 6"          
## [11] "Citrullus_lanatus_Ren2015 3"               
## [12] "Citrullus_lanatus_Ren2015 7"               
## [13] "Citrullus_lanatus_Ren2015 10"              
## [14] "Citrus_sinensis_Huang2018 4"               
## [15] "Coffea_canephora_Crouzillat2020 A"         
## [16] "Coffea_canephora_Crouzillat2020 K"         
## [17] "Cucumis_melo_Pereira2018 4"                
## [18] "Cucumis_melo_Pereira2018 7"                
## [19] "Cucumis_melo_Pereira2018 9"                
## [20] "Cucumis_sativus_Zhu2016 6"                 
## [21] "Cucumis_sativus_Zhu2016 5"                 
## [22] "Elaeis_guineensis_Yaakub2020 5"            
## [23] "Elaeis_guineensis_Yaakub2020 1"            
## [24] "Elaeis_guineensis_Yaakub2020 12"           
## [25] "Eucalyptus_grandis_Bertholome2015 1"       
## [26] "Eucalyptus_grandis_Bertholome2015 6"       
## [27] "Eucalyptus_grandis_Bertholome2015 7"       
## [28] "Eucalyptus_grandis_Bertholome2015 8"       
## [29] "Gossypium_hirsutum_Zhang2019 A11"          
## [30] "Lupinus_angustifolius_Zhou2017 8"          
## [31] "Lupinus_angustifolius_Zhou2017 10"         
## [32] "Lupinus_angustifolius_Zhou2017 19"         
## [33] "Lupinus_angustifolius_Zhou2017 15"         
## [34] "Lupinus_angustifolius_Zhou2017 3"          
## [35] "Lupinus_angustifolius_Zhou2017 7"          
## [36] "Malus_domestica_DiPierro2016 3"            
## [37] "Manihot_esculenta_ICGMC2015 4"             
## [38] "Manihot_esculenta_ICGMC2015 6"             
## [39] "Manihot_esculenta_ICGMC2015 8"             
## [40] "Manihot_esculenta_ICGMC2015 12"            
## [41] "Momordica_charantia_Mastumura2020 2"       
## [42] "Nelumbo_nucifera_Gui2018 1"                
## [43] "Nelumbo_nucifera_Gui2018 5"                
## [44] "Oryza_nivara_Ma2016 9"                     
## [45] "Oryza_nivara_Ma2016 12"                    
## [46] "Oryza_sativa_DeLeon2016 3"                 
## [47] "Oryza_sativa_DeLeon2016 5"                 
## [48] "Oryza_sativa_DeLeon2016 6"                 
## [49] "Phaseolus_vulgaris_Song2015 4"             
## [50] "Phaseolus_vulgaris_Song2015 10"            
## [51] "Phaseolus_vulgaris_Song2015 11"            
## [52] "Prunus_persica_Verde2017 3"                
## [53] "Sesamum_indicum_Wang2016 1"                
## [54] "Sesamum_indicum_Wang2016 3"                
## [55] "Sesamum_indicum_Wang2016 4"                
## [56] "Solanum_lycopersicum_Gonda2018 12"         
## [57] "Solanum_tuberosum_Endelman2016 1"          
## [58] "Solanum_tuberosum_Endelman2016 4"          
## [59] "Sorghum_bicolor_Zou2012 6"                 
## [60] "Sorghum_bicolor_Zou2012 7"                 
## [61] "Sorghum_bicolor_Zou2012 10"                
## [62] "Theobroma_cacao_Royaert2016 2"             
## [63] "Theobroma_cacao_Royaert2016 6"             
## [64] "Theobroma_cacao_Royaert2016 7"             
## [65] "Triticum_aestivum_GutierrezGonzalez2019 4B"
## [66] "Triticum_aestivum_GutierrezGonzalez2019 4D"
## [67] "Vitis_vinifera_Brault2020 4"               
## [68] "Vitis_vinifera_Brault2020 6"               
## [69] "Vitis_vinifera_Brault2020 9"               
## [70] "Vitis_vinifera_Brault2020 12"
```

```r
# The distribution of differences Values below 0 indicates that recombination
# rates are lower on the opposite position
hist(df$false_CI_rate - df$true_CI_rate, main = "", xlab = "Differences inferred/alternative centromere position", 
    breaks = 30)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

```r
# How many of wrongly oriented CI have a small difference between true and false
# position? e.g. less than 0.5 cM/Mb
threshold = 1
sum(df$true_CI_rate > df$false_CI_rate & (df$false_CI_rate - df$true_CI_rate) > -threshold, 
    na.rm = TRUE)
```

```
## [1] 45
```

```r
sum(df$true_CI_rate > df$false_CI_rate & (df$false_CI_rate - df$true_CI_rate) > -threshold, 
    na.rm = TRUE)/sum(df$true_CI_rate > df$false_CI_rate, na.rm = TRUE)
```

```
## [1] 0.6428571
```



# Smaller chromosomes recombine more than larger ones

Test the correlation between chromosome size and recombination rate.


```r
# Sample size
nrow(chromosome.stats)
```

```
## [1] 665
```


How many chromosomes with less than one CO?


```r
length(which(chromosome.stats$linkage.map.length < 50))
```

```
## [1] 11
```

```r
# proportion
length(which(chromosome.stats$linkage.map.length < 50))/nrow(chromosome.stats)
```

```
## [1] 0.01654135
```

How often do chromosome arms have multiple crossovers, versus a single one?

```r
length(which(chromosome.stats$linkage.map.length > 50 & chromosome.stats$linkage.map.length < 
    100))
```

```
## [1] 234
```

```r
length(which(chromosome.stats$linkage.map.length > 100))
```

```
## [1] 419
```

```r
# proportion
length(which(chromosome.stats$linkage.map.length < 50))/nrow(chromosome.stats)
```

```
## [1] 0.01654135
```





```r
# Spearman correlation
cor.test(chromosome.stats$phys.map.length, chromosome.stats$mean.recrate, method = "spearman")
```

```
## Warning in cor.test.default(chromosome.stats$phys.map.length,
## chromosome.stats$mean.recrate, : Cannot compute exact p-value with ties
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  chromosome.stats$phys.map.length and chromosome.stats$mean.recrate
## S = 90278103, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.8419156
```

```r
# Spearman correlation with Phylogenetically Independent Contrasts
# cor_phylo(log10(mean.recrate) ~ log10(phys.map.length), data =
# phylogenetic_traits, species = phylogenetic_traits$Species, phy = tree)

# Test the correlation with PGLMM
pglmm_mod = pglmm(log10(mean.recrate) ~ log10(phys.map.length) + (1 | Species__), 
    data = phylogenetic_traits, family = "gaussian", cov_ranef = list(Species = tree))

summary(pglmm_mod)
```

```
## Linear mixed model fit by restricted maximum likelihood
## 
## Call:log10(mean.recrate) ~ log10(phys.map.length)
## <environment: 0x7f83c9ac96c8>
## 
##  logLik     AIC     BIC 
##   654.9 -1299.8 -1283.0 
## 
## Random effects:
##             Variance Std.Dev
## 1|Species   0.038892  0.1972
## 1|Species__ 0.022891  0.1513
## residual    0.005329  0.0730
## 
## Fixed effects:
##                            Value Std.Error  Zscore    Pvalue    
## (Intercept)             0.867630  0.111049   7.813 5.584e-15 ***
## log10(phys.map.length) -0.313184  0.027989 -11.190 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


### Model selection to identify the model that best fit the data.

<table>
<caption>(\#tab:unnamed-chunk-34)Model selection based on AIC/BIC</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> model </th>
   <th style="text-align:center;"> formula </th>
   <th style="text-align:center;"> AIC </th>
   <th style="text-align:center;"> BIC </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> LM </td>
   <td style="text-align:center;"> log10(mean.recrate) ~ log10(phys.map.length) </td>
   <td style="text-align:center;"> -388.852 </td>
   <td style="text-align:center;"> -375.3527 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> LMER </td>
   <td style="text-align:center;"> log10(mean.recrate) ~ log10(phys.map.length) + (1 | Species) </td>
   <td style="text-align:center;"> -1275.590 </td>
   <td style="text-align:center;"> -1257.5906 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> LMER </td>
   <td style="text-align:center;"> log10(mean.recrate) ~ log10(phys.map.length) + (log10(phys.map.length) | Species) </td>
   <td style="text-align:center;"> -1347.495 </td>
   <td style="text-align:center;"> -1320.4960 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> PGLMM </td>
   <td style="text-align:center;"> log10(mean.recrate) ~ log10(phys.map.length) + (1|Species__) </td>
   <td style="text-align:center;"> -1286.970 </td>
   <td style="text-align:center;"> -1273.5493 </td>
  </tr>
</tbody>
</table>



```r
write.xlsx(x = model_selection, file = paste(wd, "tables/article_one_supp/tables_supplementary.xls", 
    sep = ""), sheetName = "S5_ModelSelection", row.names = FALSE, append = TRUE)
```





### Representing the relationship between chromosome size and recombination rate.




```
## 
## Call:
## lm(formula = log10(chromosome.stats$mean.recrate) ~ log10(phys.map.length), 
##     data = chromosome.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.68324 -0.14576  0.01838  0.13625  0.48238 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)             1.90461    0.02717   70.09   <2e-16 ***
## log10(phys.map.length) -0.90668    0.01569  -57.80   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1801 on 663 degrees of freedom
## Multiple R-squared:  0.8344,	Adjusted R-squared:  0.8342 
## F-statistic:  3341 on 1 and 663 DF,  p-value: < 2.2e-16
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: log10(mean.recrate) ~ log10(phys.map.length) + (1 | species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: -1283.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.5275 -0.5818 -0.0144  0.5723  3.2839 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  species  (Intercept) 0.098403 0.3137  
##  Residual             0.005344 0.0731  
## Number of obs: 665, groups:  species, 57
## 
## Fixed effects:
##                         Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)              0.89480    0.06329 193.45259   14.14   <2e-16 ***
## log10(phys.map.length)  -0.32494    0.02770 662.81832  -11.73   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## lg10(phy..) -0.753
```

```
## Warning: 'r.squaredGLMM' now calculates a revised statistic. See the help page.
```

```
##            R2m       R2c
## [1,] 0.1680834 0.9571483
```

#### Figure 1


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-37-1.png" alt="(ref:meanrecrate)"  />
<p class="caption">(\#fig:unnamed-chunk-37)(ref:meanrecrate)</p>
</div>

(ref:meanrecrate) (ref:meanrecrate) Negative correlation of recombination rates (cM/Mb, log scale) with chromosome size (Mb, log scale). Recombination rates were estimated with loess regression in windows of 100kb. Each point represents a chromosome (n=665). Species are presented in different colors (57 species). (A) The black solid line represents the linear model regression line fitted to the data. The expected simulated regression line is in red, assuming on average two crossovers per chromosome for the solid line (100cM), within an interval of 1 to 3 crossovers per chromosome (respectively lower and upper dashed red lines). (B) Within species relationships between recombination rates and chromosome size. The black dashed line represents the selected Linear Mixed Model (LMER). Colored lines for species random regressions.





### Relationship between linkage map length and chromosome length. Linkage map length is a measure of the absolute number of crossover per chromosome.


```r
cor.test(chromosomes$phys.map.length, chromosomes$linkage.map.length.correctedHW, 
    method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  chromosomes$phys.map.length and chromosomes$linkage.map.length.correctedHW
## S = 40018775, p-value = 1.899e-06
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.1835096
```









```
## 
## Call:
## lm(formula = linkage.map.length.correctedHW ~ phys.map.length, 
##     data = chromosome.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -88.612 -37.943  -7.066  28.798 236.029 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     116.67034    2.07605   56.20  < 2e-16 ***
## phys.map.length   0.05707    0.01121    5.09 4.67e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 46.04 on 663 degrees of freedom
## Multiple R-squared:  0.03761,	Adjusted R-squared:  0.03616 
## F-statistic: 25.91 on 1 and 663 DF,  p-value: 4.665e-07
```

```
## [1] 0.03761008
## attr(,"adj.r.squared")
## [1] 0.03761108
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: linkage.map.length.correctedHW ~ phys.map.length + (phys.map.length |  
##     species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 6226.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.3021 -0.5737 -0.0653  0.4890  4.9865 
## 
## Random effects:
##  Groups   Name            Variance Std.Dev. Corr 
##  species  (Intercept)     1667.311 40.833        
##           phys.map.length    3.557  1.886   -0.39
##  Residual                  389.668 19.740        
## Number of obs: 665, groups:  species, 57
## 
## Fixed effects:
##                 Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)      41.8231     6.8216 56.3289   6.131 9.12e-08 ***
## phys.map.length   2.0955     0.2736 44.4002   7.659 1.18e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## phys.mp.lng -0.489
```

```
##            R2m       R2c
## [1,] 0.4850092 0.9983054
```




```
## 
## Call:
## lm(formula = excessCO ~ chrrelativesize, data = chromosome.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -95.996 -34.787  -5.097  28.758 210.634 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -4.119      7.715  -0.534    0.594    
## chrrelativesize   76.180      7.526  10.122   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 43.68 on 663 degrees of freedom
## Multiple R-squared:  0.1338,	Adjusted R-squared:  0.1325 
## F-statistic: 102.5 on 1 and 663 DF,  p-value: < 2.2e-16
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: excessCO ~ chrrelativesize + (chrrelativesize | species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 6117.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.3118 -0.5617 -0.0818  0.4922  5.0704 
## 
## Random effects:
##  Groups   Name            Variance Std.Dev. Corr 
##  species  (Intercept)     1663.9   40.79         
##           chrrelativesize 2708.3   52.04    -0.61
##  Residual                  384.9   19.62         
## Number of obs: 665, groups:  species, 57
## 
## Fixed effects:
##                 Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)      -13.260      6.925  55.336  -1.915   0.0607 .  
## chrrelativesize   86.551      8.104  49.811  10.680 1.76e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## chrrelatvsz -0.726
```

```
##            R2m       R2c
## [1,] 0.1411685 0.8570205
```

#### Figure 2


```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: excessCO ~ chrrelativesize + (chrrelativesize | species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 6117.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.3118 -0.5617 -0.0818  0.4922  5.0704 
## 
## Random effects:
##  Groups   Name            Variance Std.Dev. Corr 
##  species  (Intercept)     1663.9   40.79         
##           chrrelativesize 2708.3   52.04    -0.61
##  Residual                  384.9   19.62         
## Number of obs: 665, groups:  species, 57
## 
## Fixed effects:
##                 Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)      -13.260      6.925  55.336  -1.915   0.0607 .  
## chrrelativesize   86.551      8.104  49.811  10.680 1.76e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## chrrelatvsz -0.726
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-45-1.png" alt="(ref:linkagemaplength-chrsize)"  />
<p class="caption">(\#fig:unnamed-chunk-45-1)(ref:linkagemaplength-chrsize)</p>
</div><div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-45-2.png" alt="(ref:linkagemaplength-chrsize)"  />
<p class="caption">(\#fig:unnamed-chunk-45-2)(ref:linkagemaplength-chrsize)</p>
</div>

(ref:linkagemaplength-chrsize) (ref:linkagemaplength-chrsize) Relationship between linkage map length and chromosome length. (A) Correlation between chromosome length (Mb) and linkage map length (cM). Each point represents a chromosome (n=665). For clarity, we substracted 50cM to the linkage map length to remove the mandatory crossover and represent only the excess of crossovers. Species are presented in different colors (57 species). The linear regression is the solid black line. The fixed regression of the Linear Mixed Model is the dashed black line. Species random slopes are in colors. Isolines of the Genome-wide Recombination Rate (GwRR) were plotted in dotted red lines to represent regions of equal recombination rates. (B) Species random intercepts as a function of the specific mean chromosome size (Mb). (C) Species random slopes as a function of the specific mean chromosome size (Mb). (D) Within species effect of relative chromosome size on linkage map length. Each point represents a chromosome (n=665). Species are presented in different colors (57 species). Black solid line is the linear regression across species while dashed black line is the linear mixed regression with a species random effect. Colored solid lines represent the Linear Mixed Model regression line (LMER) fitted to the data, with random regressions fitted to species.





#### Figure S4






Yet, it is difficult to imagine a biological mechanism that is able to consider relative sizes instead of absolute sizes. Given that number of genes are roughly constant across species and assuming that genes are distributed proportionnally to relative chromosome sizes within species, we can imagine that number of genes of a chromosome could be a good proxy of the relative size of the chromosome.

Yet, number of genes is not constant across species in our dataset, as expected if annotations are not standard and if some overlapping genes/alternative variants increases the number of gene copies. Besides, we have unequal sampling across species and for some species, we have missing chromosomes.

Hence we decided to divide gene count by the mean gene count of the species (i.e. a relative gene count) as a good proxy of relative chromosome size and compared with standardizing by the total gene count (despite uneven sampling.

Finally we kept the standardization by the total gene count.


<table>
<caption>(\#tab:unnamed-chunk-48)Gene count per species</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> Species </th>
   <th style="text-align:center;"> Number of chromosomes </th>
   <th style="text-align:center;"> Mean gene count </th>
   <th style="text-align:center;"> Total number of genes </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 6294 </td>
   <td style="text-align:center;"> 12587 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 5489 </td>
   <td style="text-align:center;"> 27445 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 4228 </td>
   <td style="text-align:center;"> 33823 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 3592 </td>
   <td style="text-align:center;"> 64656 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 5585 </td>
   <td style="text-align:center;"> 16756 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 5398 </td>
   <td style="text-align:center;"> 64778 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 4845 </td>
   <td style="text-align:center;"> 33915 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 4454 </td>
   <td style="text-align:center;"> 66816 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 2261 </td>
   <td style="text-align:center;"> 11306 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 3351 </td>
   <td style="text-align:center;"> 23457 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 2316 </td>
   <td style="text-align:center;"> 25478 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 3397 </td>
   <td style="text-align:center;"> 23780 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 1532 </td>
   <td style="text-align:center;"> 29116 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 1621 </td>
   <td style="text-align:center;"> 27554 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 976 </td>
   <td style="text-align:center;"> 13657 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 1610 </td>
   <td style="text-align:center;"> 17706 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 3102 </td>
   <td style="text-align:center;"> 34121 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 2715 </td>
   <td style="text-align:center;"> 51577 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> 26 </td>
   <td style="text-align:center;"> 2650 </td>
   <td style="text-align:center;"> 68900 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 2914 </td>
   <td style="text-align:center;"> 37885 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 2973 </td>
   <td style="text-align:center;"> 38648 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 5368 </td>
   <td style="text-align:center;"> 37577 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> 25 </td>
   <td style="text-align:center;"> 1614 </td>
   <td style="text-align:center;"> 40359 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 1475 </td>
   <td style="text-align:center;"> 29495 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 2485 </td>
   <td style="text-align:center;"> 42250 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 1772 </td>
   <td style="text-align:center;"> 31901 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 3026 </td>
   <td style="text-align:center;"> 36313 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 3154 </td>
   <td style="text-align:center;"> 37849 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 3742 </td>
   <td style="text-align:center;"> 33680 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 2545 </td>
   <td style="text-align:center;"> 27996 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 2784 </td>
   <td style="text-align:center;"> 19487 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 3351 </td>
   <td style="text-align:center;"> 26811 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 1707 </td>
   <td style="text-align:center;"> 22193 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 3807 </td>
   <td style="text-align:center;"> 34264 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 2827 </td>
   <td style="text-align:center;"> 33925 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 3124 </td>
   <td style="text-align:center;"> 37482 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 3403 </td>
   <td style="text-align:center;"> 34027 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 2899 </td>
   <td style="text-align:center;"> 28993 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 4989 </td>
   <td style="text-align:center;"> 99771 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 2863 </td>
   <td style="text-align:center;"> 31496 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 1399 </td>
   <td style="text-align:center;"> 26588 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 3900 </td>
   <td style="text-align:center;"> 39005 </td>
  </tr>
</tbody>
</table>






```
## 
## Call:
## lm(formula = excessCO ~ genecount_std, data = chromosome.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -90.982 -39.667  -0.738  28.957 152.449 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     55.453      3.679  15.072  < 2e-16 ***
## genecount_std  217.282     38.156   5.695 2.09e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 44.55 on 513 degrees of freedom
##   (150 observations deleted due to missingness)
## Multiple R-squared:  0.05945,	Adjusted R-squared:  0.05762 
## F-statistic: 32.43 on 1 and 513 DF,  p-value: 2.085e-08
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: excessCO ~ genecount_std + (1 | species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 4682.3
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.6265 -0.5573 -0.0369  0.4918  4.6119 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  species  (Intercept) 3917.9   62.59   
##  Residual              364.4   19.09   
## Number of obs: 515, groups:  species, 42
## 
## Fixed effects:
##               Estimate Std. Error      df t value Pr(>|t|)    
## (Intercept)     -9.744     10.428  46.530  -0.934    0.355    
## genecount_std  789.149     35.573 512.872  22.184   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## genecnt_std -0.366
```

```
##            R2m      R2c
## [1,] 0.2783178 0.938588
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-51-1.png" alt="(ref:linkagemaplength-genecount)"  />
<p class="caption">(\#fig:unnamed-chunk-51)(ref:linkagemaplength-genecount)</p>
</div>

(ref:linkagemaplength-genecount) (ref:linkagemaplength-genecount) Relationship between linkage map length and relative number of genes. Each point represents a chromosome (n=515). Species are presented in different colors (42 species). Black solid line is the linear regression across species while dashed black line is the linear mixed regression with a species random effect. Colored solid lines represent the Linear Mixed Model regression line (LMER) fitted to the data, with random regressions fitted to species.





## A model explaining variation in linkage map length by chromosome relative size

<!-- Genetic-map_{i,j} - 50 = b_i x Physi_map_{i,j} / GenomeSize_i -->
<!-- où b_i est le facteur d’interférence espèce spécifique. -->
<!-- mydata$RelChromSize <- mydata$phys.map.length/mydata$GenomeSize -->
<!-- lmer( (mydata$linkage.map.length - 50) ~ 0 + mydata$RelChromSize|mydata$set ) -->
<!-- Le R2 est de 0.95 -->
<!-- Si on rajoute un intercept random ça n’améliore quasiment pas -->
<!-- lmer( (mydata$linkage.map.length-50) ~ 1|mydata$set + mydata$RelChromSize|mydata$set ) -->
<!-- R2 = 0.962 -->


Ultimately, we can think of a model:

$$ d = 50 + a*n_i $$
$$ d = 50 + a*n*\frac{L_i}{L} $$
Where $n$ is the number of genes of the species, $n_i$ the number of genes in chromosome $i$, $L$ the total physical length of the genome and $L_i$ the physical length of chromosome $i$, and, finally, $a$  is a coefficient supposedly constant across species.

<!-- A species has N chromosomes of length Li (i from 1 to N). We assume that the number of genes is proportional to gene length: -->
<!-- So ni = n x Li / L where L is the total genome size and n is the toal number of genes. -->
<!-- If we assume that the excess of genetic distance, beyond the obligate 50 cM, is proportional to the number of genes, not the physical length, the genetic distance is given by: -->
<!-- di = 50 + a x ni -->
<!-- di = 50 + a x n x Li / L  -->
<!-- If we assume that the number of genes, n, and the a coefficient are constant across species we have the following schematic predictions: -->







# Intra-chromosomal heterogeneity in recombination

We observed contrasted patterns between species.

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-52-1.png" alt="(ref:heterogeneity-recombination-maps)"  />
<p class="caption">(\#fig:unnamed-chunk-52)(ref:heterogeneity-recombination-maps)</p>
</div>

(ref:heterogeneity-recombination-maps) (ref:heterogeneity-recombination-maps) The diversity of recombination landscapes in angiosperms demonstrated in six different species. Recombination landscapes are similar within species (pooled chromosomes). Physical distances were scaled for comparison of chromosomes with different sizes. Estimates of recombination rates were obtained by 1,000 bootstraps of loci in windows of 100kb with loess regression and automatic span calibration. The mean species landscape (dashed line) was estimated by computing the mean recombination rate in 100 bins along the chromosome axis, all chromosomes pooled. Similarly, the lower and upper boundaries (pale ribbon) were estimated by taking the minimum and maximum recombination rates in 100 bins. One chromosome is represented in a solid line for each species, with the physical position of the centromere pointed by a dot. Species ordered by ascending mean chromosome size. (A) *Capsella rubella* chromosome 1 and n = 7 chromosomes pooled (B) *Arabidopsis thaliana* chromosome 1 and n = 5 chromosomes pooled (C) *Malus domestica* chromosome 2 and n = 17 chromosomes pooled (D) *Eucalyptus grandis* chromosome 10  and n = 11 chromosomes pooled (E) *Nelumbo nucifera* chromosome 3 and n = 8 chromosomes pooled (F) *Zea mays* chromosome 10 and n = 10 chromosomes pooled.


#### Figure 3









## Broken stick model

We used the broken stick model (k = 10) to visually describe the diversity of patterns.











<!-- (ref:broken-stick-annotated) (ref:broken-stick-annotated) Structure of the recombination within chromosomes. Relative recombination rates along the chromosome estimated in ten bins under the broken stick model, with species ordered by ascending variance of relative recombination rates (n=665). (A) The relative recombination rate is the log-transformed ratio of the expected relative genomic length (one tenth) divided by observed relative genomic length. It means that values below zero are recombination rates lower than expected while values above zero are recombination rates higher than expected. When information is available, the centromere position on the chromosome is mapped as a red dot and chromosomes are oriented with the longer arm on top. (B) Chromosome length measured in Mb. (C) Genome-wide recombination rate, i.e. the mean recombination rate estimated in 100kb windows (cM/Mb). -->



#### Figure 4

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-58-1.png" alt="(ref:broken-stick)"  />
<p class="caption">(\#fig:unnamed-chunk-58-1)(ref:broken-stick)</p>
</div><div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-58-2.png" alt="(ref:broken-stick)"  />
<p class="caption">(\#fig:unnamed-chunk-58-2)(ref:broken-stick)</p>
</div>

(ref:broken-stick) (ref:broken-stick) Structure of the recombination within chromosomes. relative recombination rates along the chromosome estimated in ten bins under the broken stick model, with species ordered by ascending variance of relative recombination rates (n=665). The relative recombination rate is the log-transformed ratio of the expected relative genetic length (one tenth) divided by observed relative genetic length. It means that values below zero are recombination rates lower than expected while values above zero are recombination rates higher than expected. When information is available, the centromere position on the chromosome is mapped as a red dot and chromosomes are oriented with the longer arm on top.





## Heterogeneity measures

### Comparative approach - choice of the statistic

Three statistics for heterogeneity:
* coefficient of variation
* Gini coefficient
* variance in broken stick proportions

Model comparison approach:


```r
lm_model = lm(cv.recrate ~ phys.map.length, data = chromosome.stats)
AIC(lm_model)
```

```
## [1] 754.5539
```

```r
BIC(lm_model)
```

```
## [1] 768.0533
```

```r
r.squaredLR(lm_model)
```

```
## [1] 0.1728832
## attr(,"adj.r.squared")
## [1] 0.2362929
```

```r
# par(mfrow = c(2,2)) plot(lm_model, which = c(1,2,3,5)) par(mfrow = c(1,1))

lm_model = lm(gini ~ phys.map.length, data = chromosome.stats)
AIC(lm_model)
```

```
## [1] -637.5537
```

```r
BIC(lm_model)
```

```
## [1] -624.0544
```

```r
r.squaredLR(lm_model)
```

```
## [1] 0.1801051
## attr(,"adj.r.squared")
## [1] -0.1555344
```

```r
# par(mfrow = c(2,2)) plot(lm_model, which = c(1,2,3,5)) par(mfrow = c(1,1))

lm_model = lm(brokenstick_pvariance_scaled ~ phys.map.length, data = chromosome.stats)
AIC(lm_model)
```

```
## [1] 1758.335
```

```r
BIC(lm_model)
```

```
## [1] 1771.775
```

```r
r.squaredLR(lm_model)
```

```
## [1] 0.1381831
## attr(,"adj.r.squared")
## [1] 0.1467909
```

```r
# par(mfrow = c(2,2)) plot(lm_model, which = c(1,2,3,5)) par(mfrow = c(1,1))
```


```r
lm_model = lm(cv.recrate ~ phys.map.length, data = species.stats)
AIC(lm_model)
```

```
## [1] 57.17865
```

```r
BIC(lm_model)
```

```
## [1] 63.3078
```

```r
r.squaredLR(lm_model)
```

```
## [1] 0.2443293
## attr(,"adj.r.squared")
## [1] 0.3530197
```

```r
lm_model = lm(pvariance ~ phys.map.length, data = species.stats)
AIC(lm_model)
```

```
## [1] 142.3715
```

```r
BIC(lm_model)
```

```
## [1] 148.5007
```

```r
r.squaredLR(lm_model)
```

```
## [1] 0.2263934
## attr(,"adj.r.squared")
## [1] 0.2436196
```



### Larger chromosomes have higher heterogeneity?

To test the hypothesis that larger chromosomes have a higher heterogeneity, we used the variance in broken stick segment sizes as a proxy of recombination heterogeneity. But we need to correct for species random effect.


### Broken stick variance


```r
# Distribution of variables
hist(chromosome.stats$brokenstick_pvariance_scaled, breaks = 30, main = "", xlab = "Scaled broken stick variance")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

Across species: Chromosomes pooled per species. Is there a common pattern across species?


```r
nrow(species.stats)
```

```
## [1] 57
```

```r
# Model selection Linear regression
lm_model = lm(pvariance ~ phys.map.length, data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = pvariance ~ phys.map.length, data = species.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.91315 -0.53797 -0.01841  0.60257  1.73917 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     -0.3064077  0.1279780  -2.394 0.020092 *  
## phys.map.length  0.0025801  0.0006431   4.012 0.000183 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8148 on 55 degrees of freedom
## Multiple R-squared:  0.2264,	Adjusted R-squared:  0.2123 
## F-statistic:  16.1 on 1 and 55 DF,  p-value: 0.0001835
```

```r
AIC(lm_model)
```

```
## [1] 142.3715
```

```r
BIC(lm_model)
```

```
## [1] 148.5007
```

```r
# Quadratic regression
lm_model = lm(pvariance ~ phys.map.length + I(phys.map.length^2), data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = pvariance ~ phys.map.length + I(phys.map.length^2), 
##     data = species.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -2.05672 -0.42745  0.00285  0.56475  1.65223 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)   
## (Intercept)          -5.621e-01  1.746e-01  -3.219  0.00218 **
## phys.map.length       8.273e-03  2.803e-03   2.952  0.00466 **
## I(phys.map.length^2) -8.921e-06  4.281e-06  -2.084  0.04193 * 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.7911 on 54 degrees of freedom
## Multiple R-squared:  0.284,	Adjusted R-squared:  0.2574 
## F-statistic: 10.71 on 2 and 54 DF,  p-value: 0.0001211
```

```r
AIC(lm_model)
```

```
## [1] 139.9634
```

```r
BIC(lm_model)
```

```
## [1] 148.1356
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-63-1.png" alt="(ref:brokenstickvariance-species)"  />
<p class="caption">(\#fig:unnamed-chunk-63)(ref:brokenstickvariance-species)</p>
</div>
(ref:brokenstickvariance-species) The broken stick variance depends on the mean chromosome size across species (n = 38). The broken stick variance was log-transformed and scaled for comparison between species.


Within species: Chromosome-level correlation.


```
## 
## Call:
## lm(formula = brokenstick_pvariance_scaled ~ log10(phys.map.length), 
##     data = chromosome.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -3.08746 -0.53587  0.06989  0.62189  1.81720 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            -1.79180    0.13431  -13.34   <2e-16 ***
## log10(phys.map.length)  1.06713    0.07731   13.80   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.8801 on 650 degrees of freedom
##   (13 observations deleted due to missingness)
## Multiple R-squared:  0.2267,	Adjusted R-squared:  0.2255 
## F-statistic: 190.5 on 1 and 650 DF,  p-value: < 2.2e-16
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: 
## brokenstick_pvariance_scaled ~ log10(phys.map.length) + (log10(phys.map.length) |  
##     species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 1259
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.7931 -0.5218  0.0826  0.5785  2.9358 
## 
## Random effects:
##  Groups   Name                   Variance Std.Dev. Corr 
##  species  (Intercept)            4.4676   2.1137        
##           log10(phys.map.length) 1.8818   1.3718   -0.93
##  Residual                        0.2859   0.5347        
## Number of obs: 652, groups:  species, 57
## 
## Fixed effects:
##                         Estimate Std. Error        df t value Pr(>|t|)
## (Intercept)            -0.004562   0.426866 38.830563  -0.011    0.992
## log10(phys.map.length) -0.118446   0.266835 37.912192  -0.444    0.660
## 
## Correlation of Fixed Effects:
##             (Intr)
## lg10(phy..) -0.961
```

```
##              R2m       R2c
## [1,] 0.002007696 0.7944873
```

Indeed, the adjusted R-squared of the linear model is low and the mixed model explains a low proportion of inter-species effect (R2m, marginal $R_2$), with a non-significant coefficient.


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-65-1.png" alt="(ref:relativechromosomesize-linkagemaplength)"  />
<p class="caption">(\#fig:unnamed-chunk-65)(ref:relativechromosomesize-linkagemaplength)</p>
</div>



### Coefficient of variation


```r
# Distribution of variables
hist(chromosome.stats$cv.recrate, breaks = 30, main = "", xlab = "Coefficient of variation")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-66-1.png)<!-- -->

Across species: Chromosomes pooled per species.



```r
nrow(species.stats)
```

```
## [1] 57
```

```r
# Model selection Linear regression
lm_model = lm(cv.recrate ~ phys.map.length, data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = cv.recrate ~ phys.map.length, data = species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.6389 -0.2502 -0.1124  0.2086  1.6000 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     0.7934943  0.0606156  13.091  < 2e-16 ***
## phys.map.length 0.0012845  0.0003046   4.217 9.31e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3859 on 55 degrees of freedom
## Multiple R-squared:  0.2443,	Adjusted R-squared:  0.2306 
## F-statistic: 17.78 on 1 and 55 DF,  p-value: 9.308e-05
```

```r
AIC(lm_model)
```

```
## [1] 57.17865
```

```r
BIC(lm_model)
```

```
## [1] 63.3078
```

```r
# Quadratic regression
lm_model = lm(cv.recrate ~ phys.map.length + I(phys.map.length^2), data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = cv.recrate ~ phys.map.length + I(phys.map.length^2), 
##     data = species.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.06646 -0.20389 -0.04645  0.23105  1.17056 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           6.013e-01  7.749e-02   7.759 2.41e-10 ***
## phys.map.length       5.565e-03  1.244e-03   4.475 3.99e-05 ***
## I(phys.map.length^2) -6.707e-06  1.900e-06  -3.531 0.000857 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3511 on 54 degrees of freedom
## Multiple R-squared:  0.386,	Adjusted R-squared:  0.3633 
## F-statistic: 16.98 on 2 and 54 DF,  p-value: 1.904e-06
```

```r
AIC(lm_model)
```

```
## [1] 47.34021
```

```r
BIC(lm_model)
```

```
## [1] 55.51242
```



```
## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties

## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties
```

```
## [1] -0.1833778
```

```
##       2.5% 
## -0.2927582
```

```
##      97.5% 
## -0.0716899
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-68-1.png)<!-- -->



Selected model for simplicity


```r
# The mean species correlation estimated by 1,000 bootstraps
mean(boot)
```

```
## [1] -0.1833778
```

```r
quantile(boot, 0.025)
```

```
##       2.5% 
## -0.2927582
```

```r
quantile(boot, 0.975)
```

```
##      97.5% 
## -0.0716899
```

```r
# Species correlations

summary(species.stats$correlation_cvrecrate)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -1.00000 -0.50000 -0.25000 -0.18266  0.07143  1.00000
```

```r
cor.test(species.stats$cv.recrate, species.stats$phys.map.length, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  species.stats$cv.recrate and species.stats$phys.map.length
## S = 13928, p-value = 1.398e-05
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.5486129
```

```r
# Linear regression
lm_model = lm(cv.recrate ~ phys.map.length, data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = cv.recrate ~ phys.map.length, data = species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.6389 -0.2502 -0.1124  0.2086  1.6000 
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     0.7934943  0.0606156  13.091  < 2e-16 ***
## phys.map.length 0.0012845  0.0003046   4.217 9.31e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3859 on 55 degrees of freedom
## Multiple R-squared:  0.2443,	Adjusted R-squared:  0.2306 
## F-statistic: 17.78 on 1 and 55 DF,  p-value: 9.308e-05
```



<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-70-1.png" alt="(ref:cvrecrate-species)"  />
<p class="caption">(\#fig:unnamed-chunk-70)(ref:cvrecrate-species)</p>
</div>

(ref:cvrecrate-species) The coefficient of variation depends on the mean chromosome size across species (n = 38). The broken stick variance was log-transformed and scaled for comparison between species. The inset presents the distribution of Spearman's correlation within species between the coefficient of variation and chromosome length. Mean correlation and 95% C.I. estimated by 1,000 bootstraps.


Within species: Chromosome-level correlation.


```
## 
## Call:
## lm(formula = cv.recrate ~ log10(phys.map.length), data = chromosome.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.0652 -0.2445 -0.0515  0.1790  4.0242 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)             0.02402    0.06051   0.397    0.692    
## log10(phys.map.length)  0.54023    0.03493  15.467   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.401 on 663 degrees of freedom
## Multiple R-squared:  0.2652,	Adjusted R-squared:  0.264 
## F-statistic: 239.2 on 1 and 663 DF,  p-value: < 2.2e-16
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: cv.recrate ~ log10(phys.map.length) + (log10(phys.map.length) |  
##     species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 216.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.9743 -0.4893 -0.0659  0.4373  9.7273 
## 
## Random effects:
##  Groups   Name                   Variance Std.Dev. Corr 
##  species  (Intercept)            0.74952  0.8657        
##           log10(phys.map.length) 0.36995  0.6082   -0.96
##  Residual                        0.06028  0.2455        
## Number of obs: 665, groups:  species, 57
## 
## Fixed effects:
##                        Estimate Std. Error      df t value Pr(>|t|)  
## (Intercept)              0.4867     0.1744 24.1775   2.790   0.0101 *
## log10(phys.map.length)   0.2099     0.1145 33.0175   1.833   0.0758 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## lg10(phy..) -0.975
```

```
##             R2m       R2c
## [1,] 0.03775708 0.7398221
```

For the coefficient of variation, the adjusted R-squared of the linear model is low and the mixed model explains a low proportion of inter-species effect (R2m, marginal $R_2$), as well as for the broken stick model.


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-72-1.png" alt="(ref:chromosomesize-cvrecrate)"  />
<p class="caption">(\#fig:unnamed-chunk-72)(ref:chromosomesize-cvrecrate)</p>
</div>

(ref:chromosomesize-cvrecrate) The coefficient of variation within species does not always vary with chromosome size.




## Species-specific heterogeneity of recombination landscapes 

Test explicitly the hypothesis that recombination landscapes are species specific. The variance in broken stick proportions and coefficient of variation should be correlated to the species.


```r
# Test the species effect
aov_mod = aov(brokenstick_pvariance_scaled ~ species, data = chromosome.stats)
summary(aov_mod)
```

```
##              Df Sum Sq Mean Sq F value Pr(>F)    
## species      56  466.9   8.338   26.95 <2e-16 ***
## Residuals   595  184.1   0.309                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 13 observations deleted due to missingness
```

```r
lm_mod = lm(brokenstick_pvariance_scaled ~ species, data = chromosome.stats)
summary(lm_mod)
```

```
## 
## Call:
## lm(formula = brokenstick_pvariance_scaled ~ species, data = chromosome.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -2.81084 -0.27685  0.03808  0.32639  1.65400 
## 
## Coefficients:
##                                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                     1.31771    0.39333   3.350 0.000859 ***
## speciesArabidopsis thaliana    -2.54768    0.46539  -5.474 6.48e-08 ***
## speciesArachis duranensis       0.09747    0.44599   0.219 0.827077    
## speciesArachis hypogaea        -1.55533    0.41460  -3.751 0.000193 ***
## speciesBoechera stricta        -2.68769    0.44599  -6.026 2.94e-09 ***
## speciesBrachypodium distachyon -1.37641    0.50778  -2.711 0.006909 ** 
## speciesBrassica napus          -1.46820    0.42759  -3.434 0.000637 ***
## speciesBrassica rapa           -2.08077    0.44599  -4.665 3.81e-06 ***
## speciesCamelina sativa         -2.00421    0.42048  -4.766 2.36e-06 ***
## speciesCamellia sinensis       -2.53904    0.46539  -5.456 7.16e-08 ***
## speciesCapsella rubella        -2.17034    0.44599  -4.866 1.46e-06 ***
## speciesCapsicum annuum         -0.43555    0.42759  -1.019 0.308797    
## speciesCenchrus americanus     -0.09636    0.44599  -0.216 0.829012    
## speciesCitrullus lanatus       -1.59629    0.43087  -3.705 0.000231 ***
## speciesCitrus sinensis         -0.96389    0.43975  -2.192 0.028773 *  
## speciesCoffea canephora        -0.43031    0.42759  -1.006 0.314649    
## speciesCucumis melo            -0.73225    0.42759  -1.713 0.087323 .  
## speciesCucumis sativus         -3.46655    0.44599  -7.773 3.39e-14 ***
## speciesCucurbita maxima        -2.59422    0.41460  -6.257 7.50e-10 ***
## speciesCucurbita pepo          -2.39334    0.42048  -5.692 1.97e-08 ***
## speciesDioscorea alata         -1.03498    0.42250  -2.450 0.014586 *  
## speciesDraba nivalis           -2.03537    0.45417  -4.481 8.89e-06 ***
## speciesElaeis guineensis       -2.41353    0.43087  -5.602 3.25e-08 ***
## speciesEucalyptus grandis      -1.97041    0.42759  -4.608 4.97e-06 ***
## speciesGlycine max             -0.52113    0.41351  -1.260 0.208071    
## speciesGossypium hirsutum      -0.59942    0.40817  -1.469 0.142484    
## speciesGossypium raimondii     -0.86785    0.42250  -2.054 0.040403 *  
## speciesHelianthus annuus       -0.09566    0.42250  -0.226 0.820964    
## speciesHordeum vulgare          0.20261    0.44599   0.454 0.649786    
## speciesJuglans regia           -1.79731    0.42048  -4.274 2.23e-05 ***
## speciesLupinus albus           -1.43969    0.40876  -3.522 0.000461 ***
## speciesLupinus angustifolius   -0.52135    0.41252  -1.264 0.206792    
## speciesMalus domestica         -1.96454    0.41582  -4.724 2.88e-06 ***
## speciesMangifera indica        -1.73166    0.41252  -4.198 3.11e-05 ***
## speciesManihot esculenta       -1.98331    0.41460  -4.784 2.17e-06 ***
## speciesMomordica charantia     -1.29220    0.42759  -3.022 0.002619 ** 
## speciesNelumbo nucifera        -2.60544    0.45417  -5.737 1.54e-08 ***
## speciesOryza nivara            -1.76165    0.42484  -4.147 3.87e-05 ***
## speciesOryza sativa            -2.10905    0.42484  -4.964 9.02e-07 ***
## speciesPanicum hallii          -0.11654    0.43484  -0.268 0.788788    
## speciesPhaseolus vulgaris      -0.26853    0.42759  -0.628 0.530236    
## speciesPrunus mume             -2.56481    0.44599  -5.751 1.42e-08 ***
## speciesPrunus persica          -2.21901    0.43975  -5.046 6.00e-07 ***
## speciesQuercus sp              -2.44718    0.42484  -5.760 1.35e-08 ***
## speciesRaphanus sativus        -1.33445    0.45417  -2.938 0.003429 ** 
## speciesSesamum indicum         -2.11975    0.42250  -5.017 6.93e-07 ***
## speciesSetaria italica         -1.30421    0.43484  -2.999 0.002819 ** 
## speciesSolanum lycopersicum     0.28861    0.42484   0.679 0.497186    
## speciesSolanum tuberosum       -0.36958    0.42484  -0.870 0.384693    
## speciesSorghum bicolor         -0.66342    0.43087  -1.540 0.124155    
## speciesTheobroma cacao         -1.37108    0.43087  -3.182 0.001538 ** 
## speciesTriticum aestivum       -0.34092    0.41252  -0.826 0.408890    
## speciesTriticum dicoccoides    -0.13988    0.43087  -0.325 0.745556    
## speciesTriticum urartu         -0.60336    0.44599  -1.353 0.176614    
## speciesVigna unguiculata       -1.14472    0.42759  -2.677 0.007630 ** 
## speciesVitis vinifera          -1.69584    0.41351  -4.101 4.68e-05 ***
## speciesZea mays                -0.84711    0.43087  -1.966 0.049756 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.5562 on 595 degrees of freedom
##   (13 observations deleted due to missingness)
## Multiple R-squared:  0.7172,	Adjusted R-squared:  0.6906 
## F-statistic: 26.95 on 56 and 595 DF,  p-value: < 2.2e-16
```


```r
# Test the species effect
aov_mod = aov(cv.recrate ~ species, data = chromosome.stats)
summary(aov_mod)
```

```
##              Df Sum Sq Mean Sq F value Pr(>F)    
## species      56  105.5  1.8837   28.92 <2e-16 ***
## Residuals   608   39.6  0.0651                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
lm_mod = lm(cv.recrate ~ species, data = chromosome.stats)
summary(lm_mod)
```

```
## 
## Call:
## lm(formula = cv.recrate ~ species, data = chromosome.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.00545 -0.11579 -0.01508  0.11122  2.70661 
## 
## Coefficients:
##                                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                     1.55895    0.18047   8.638  < 2e-16 ***
## speciesArabidopsis thaliana    -1.08520    0.21353  -5.082 4.97e-07 ***
## speciesArachis duranensis       0.03361    0.20177   0.167 0.867779    
## speciesArachis hypogaea        -0.75205    0.19023  -3.953 8.61e-05 ***
## speciesBoechera stricta        -1.08469    0.20463  -5.301 1.62e-07 ***
## speciesBrachypodium distachyon -0.70161    0.23298  -3.011 0.002708 ** 
## speciesBrassica napus          -0.81326    0.19493  -4.172 3.46e-05 ***
## speciesBrassica rapa           -0.95514    0.20463  -4.668 3.75e-06 ***
## speciesCamelina sativa         -1.00457    0.19212  -5.229 2.35e-07 ***
## speciesCamellia sinensis       -1.11900    0.21353  -5.240 2.21e-07 ***
## speciesCapsella rubella        -0.90493    0.20463  -4.422 1.16e-05 ***
## speciesCapsicum annuum         -0.18093    0.19619  -0.922 0.356778    
## speciesCenchrus americanus      1.12162    0.20463   5.481 6.20e-08 ***
## speciesCitrullus lanatus       -0.83471    0.19769  -4.222 2.79e-05 ***
## speciesCitrus sinensis         -0.52343    0.20177  -2.594 0.009711 ** 
## speciesCoffea canephora        -0.39652    0.19619  -2.021 0.043710 *  
## speciesCucumis melo            -0.41467    0.19619  -2.114 0.034954 *  
## speciesCucumis sativus         -1.18774    0.20463  -5.804 1.04e-08 ***
## speciesCucurbita maxima        -1.00608    0.18973  -5.303 1.60e-07 ***
## speciesCucurbita pepo          -0.94488    0.19079  -4.952 9.51e-07 ***
## speciesDioscorea alata         -0.29643    0.19293  -1.536 0.124942    
## speciesDraba nivalis           -0.97952    0.20463  -4.787 2.13e-06 ***
## speciesElaeis guineensis       -1.00737    0.19619  -5.135 3.81e-07 ***
## speciesEucalyptus grandis      -0.85847    0.19619  -4.376 1.42e-05 ***
## speciesGlycine max             -0.43668    0.18973  -2.302 0.021695 *  
## speciesGossypium hirsutum      -0.29799    0.18728  -1.591 0.112098    
## speciesGossypium raimondii     -0.32976    0.19386  -1.701 0.089438 .  
## speciesHelianthus annuus        0.21115    0.19386   1.089 0.276499    
## speciesHordeum vulgare          0.10148    0.20463   0.496 0.620138    
## speciesJuglans regia           -0.90768    0.19212  -4.724 2.87e-06 ***
## speciesLupinus albus           -0.86176    0.18755  -4.595 5.27e-06 ***
## speciesLupinus angustifolius   -0.23931    0.18928  -1.264 0.206587    
## speciesMalus domestica         -0.99378    0.19079  -5.209 2.61e-07 ***
## speciesMangifera indica        -0.89565    0.18928  -4.732 2.77e-06 ***
## speciesManihot esculenta       -0.87401    0.19023  -4.594 5.28e-06 ***
## speciesMomordica charantia     -0.78715    0.19619  -4.012 6.77e-05 ***
## speciesNelumbo nucifera        -1.07501    0.20177  -5.328 1.40e-07 ***
## speciesOryza nivara            -0.87619    0.19493  -4.495 8.33e-06 ***
## speciesOryza sativa            -0.97567    0.19493  -5.005 7.32e-07 ***
## speciesPanicum hallii          -0.41636    0.19952  -2.087 0.037316 *  
## speciesPhaseolus vulgaris      -0.22964    0.19619  -1.171 0.242257    
## speciesPrunus mume             -1.03932    0.20463  -5.079 5.06e-07 ***
## speciesPrunus persica          -1.05811    0.20177  -5.244 2.17e-07 ***
## speciesQuercus sp              -1.01806    0.19493  -5.223 2.42e-07 ***
## speciesRaphanus sativus        -0.89503    0.20839  -4.295 2.03e-05 ***
## speciesSesamum indicum         -0.96275    0.19386  -4.966 8.88e-07 ***
## speciesSetaria italica         -0.64874    0.19952  -3.252 0.001211 ** 
## speciesSolanum lycopersicum     0.05541    0.19493   0.284 0.776311    
## speciesSolanum tuberosum       -0.50016    0.19493  -2.566 0.010531 *  
## speciesSorghum bicolor         -0.60795    0.19769  -3.075 0.002198 ** 
## speciesTheobroma cacao         -0.75106    0.19769  -3.799 0.000160 ***
## speciesTriticum aestivum       -0.09606    0.18928  -0.508 0.611966    
## speciesTriticum dicoccoides    -0.02178    0.19769  -0.110 0.912328    
## speciesTriticum urartu         -0.53515    0.20463  -2.615 0.009139 ** 
## speciesVigna unguiculata       -0.75591    0.19619  -3.853 0.000129 ***
## speciesVitis vinifera          -0.77720    0.18973  -4.096 4.77e-05 ***
## speciesZea mays                -0.40769    0.19769  -2.062 0.039611 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.2552 on 608 degrees of freedom
## Multiple R-squared:  0.727,	Adjusted R-squared:  0.7019 
## F-statistic: 28.92 on 56 and 608 DF,  p-value: < 2.2e-16
```





## COs are biased toward the periphery for many species

Estimate the bias toward the periphery.


```r
# Estimate the mean periphery-bias ratio with bootstrapped C.I.
nboot = 1000
boot = numeric(nboot)
for (i in 1:nboot) {
    boot[i] = mean(sample(chromosome.stats$peripherybias_ratio, replace = TRUE), 
        na.rm = TRUE)
}
mean(boot)
```

```
## [1] 2.187515
```

```r
quantile(boot, 0.025)
```

```
##     2.5% 
## 2.068109
```

```r
quantile(boot, 0.975)
```

```
##    97.5% 
## 2.329313
```

```r
# Results
hist(chromosome.stats$peripherybias_ratio, xlab = "Periphery-bias ratio", breaks = 30, 
    main = "")
abline(v = 1, col = "black")
abline(v = mean(boot), col = "red", lwd = 2)
abline(v = quantile(boot, 0.025), lty = "dashed", col = "red", lwd = 2)
abline(v = quantile(boot, 0.975), lty = "dashed", col = "red", lwd = 2)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-75-1.png)<!-- -->

```r
# How many maps have a periphery-bias ratio over 1 (i.e. more recombination at
# the tips than at the center)?
sum(chromosome.stats$peripherybias_ratio > 1, na.rm = TRUE)
```

```
## [1] 484
```

```r
# And in proportion
sum(chromosome.stats$peripherybias_ratio > 1, na.rm = TRUE)/nrow(chromosome.stats)
```

```
## [1] 0.7278195
```

```r
# How many maps have a periphery-bias ratio under 1 (i.e. more recombination at
# the center than on the tips)?
sum(chromosome.stats$peripherybias_ratio < 1, na.rm = TRUE)
```

```
## [1] 169
```

```r
# And in proportion
sum(chromosome.stats$peripherybias_ratio < 1, na.rm = TRUE)/nrow(chromosome.stats)
```

```
## [1] 0.2541353
```
Sample size. 


```r
# Sample size
nrow(species.stats)
```

```
## [1] 57
```

Test which species have a significant periphery-bias ratio (supplementary table).



```
## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties

## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties

## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties

## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties

## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties

## Warning in
## cor.test.default(chromosome.stats$phys.map.length[which(chromosome.stats$species
## == : Cannot compute exact p-value with ties
```

```
## [1] 0.1201095
```

```
##       2.5% 
## 0.01758947
```

```
##     97.5% 
## 0.2195359
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-77-1.png)<!-- -->



```r
kable(species.stats[, c(1:4, 9:10)], caption = "Species correlation between the periphery-bias ratio and chromosome length.", 
    labels = c("Species", "Chromosome number", "Chromosome size", "Mean periphery-bias", 
        "Spearman correlation", "p-value"), align = "c")
```

<table>
<caption>(\#tab:unnamed-chunk-78)Species correlation between the periphery-bias ratio and chromosome length.</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> species </th>
   <th style="text-align:center;"> chromosome_number </th>
   <th style="text-align:center;"> phys.map.length </th>
   <th style="text-align:center;"> peripherybias_ratio </th>
   <th style="text-align:center;"> correlation_peripherybias </th>
   <th style="text-align:center;"> correlation_peripherybias_p </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 517.44324 </td>
   <td style="text-align:center;"> 5.7567440 </td>
   <td style="text-align:center;"> -1.0000000 </td>
   <td style="text-align:center;"> 1.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 23.82927 </td>
   <td style="text-align:center;"> 1.6017099 </td>
   <td style="text-align:center;"> -0.7000000 </td>
   <td style="text-align:center;"> 0.2333333 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 100.95097 </td>
   <td style="text-align:center;"> 4.6362972 </td>
   <td style="text-align:center;"> 0.8095238 </td>
   <td style="text-align:center;"> 0.0217758 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 124.68497 </td>
   <td style="text-align:center;"> 2.0663790 </td>
   <td style="text-align:center;"> 0.1682147 </td>
   <td style="text-align:center;"> 0.5032006 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Boechera stricta </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 26.52555 </td>
   <td style="text-align:center;"> 1.5082194 </td>
   <td style="text-align:center;"> -0.0714286 </td>
   <td style="text-align:center;"> 0.9063492 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 45.62172 </td>
   <td style="text-align:center;"> 2.6055663 </td>
   <td style="text-align:center;"> 0.5000000 </td>
   <td style="text-align:center;"> 1.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 44.71564 </td>
   <td style="text-align:center;"> 1.6530123 </td>
   <td style="text-align:center;"> 0.2242424 </td>
   <td style="text-align:center;"> 0.5366881 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 26.14066 </td>
   <td style="text-align:center;"> 1.5540350 </td>
   <td style="text-align:center;"> 0.0000000 </td>
   <td style="text-align:center;"> 1.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 30.63062 </td>
   <td style="text-align:center;"> 1.3559678 </td>
   <td style="text-align:center;"> -0.1118881 </td>
   <td style="text-align:center;"> 0.7327754 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 222.10960 </td>
   <td style="text-align:center;"> 1.1966508 </td>
   <td style="text-align:center;"> 0.5000000 </td>
   <td style="text-align:center;"> 1.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 15.83906 </td>
   <td style="text-align:center;"> 0.6845333 </td>
   <td style="text-align:center;"> -0.1428571 </td>
   <td style="text-align:center;"> 0.8027778 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsicum annuum </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 227.12175 </td>
   <td style="text-align:center;"> 3.1911317 </td>
   <td style="text-align:center;"> -0.2272727 </td>
   <td style="text-align:center;"> 0.5030670 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cenchrus americanus </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 223.50536 </td>
   <td style="text-align:center;"> 5.3945568 </td>
   <td style="text-align:center;"> 0.5357143 </td>
   <td style="text-align:center;"> 0.2357143 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrullus lanatus </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 33.59728 </td>
   <td style="text-align:center;"> 1.3777785 </td>
   <td style="text-align:center;"> 0.6121212 </td>
   <td style="text-align:center;"> 0.0664691 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 26.02017 </td>
   <td style="text-align:center;"> 2.1229883 </td>
   <td style="text-align:center;"> -0.1666667 </td>
   <td style="text-align:center;"> 0.7033234 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Coffea canephora </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 53.00321 </td>
   <td style="text-align:center;"> 2.8290326 </td>
   <td style="text-align:center;"> 0.0727273 </td>
   <td style="text-align:center;"> 0.8388249 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 30.99121 </td>
   <td style="text-align:center;"> 2.6899882 </td>
   <td style="text-align:center;"> 0.4363636 </td>
   <td style="text-align:center;"> 0.1825319 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 27.40843 </td>
   <td style="text-align:center;"> 0.6934489 </td>
   <td style="text-align:center;"> 0.7142857 </td>
   <td style="text-align:center;"> 0.0880952 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 10.34886 </td>
   <td style="text-align:center;"> 0.6825170 </td>
   <td style="text-align:center;"> 0.1053922 </td>
   <td style="text-align:center;"> 0.6872983 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 10.78700 </td>
   <td style="text-align:center;"> 0.9581566 </td>
   <td style="text-align:center;"> 0.0441176 </td>
   <td style="text-align:center;"> 0.8685676 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 23.84612 </td>
   <td style="text-align:center;"> 3.2114984 </td>
   <td style="text-align:center;"> -0.1384615 </td>
   <td style="text-align:center;"> 0.6375158 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Draba nivalis </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 34.11541 </td>
   <td style="text-align:center;"> 1.4252100 </td>
   <td style="text-align:center;"> 0.1785714 </td>
   <td style="text-align:center;"> 0.7130952 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 46.94000 </td>
   <td style="text-align:center;"> 0.8058865 </td>
   <td style="text-align:center;"> 0.0181818 </td>
   <td style="text-align:center;"> 0.9728412 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 55.68998 </td>
   <td style="text-align:center;"> 1.1661899 </td>
   <td style="text-align:center;"> 0.3818182 </td>
   <td style="text-align:center;"> 0.2483806 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 48.12719 </td>
   <td style="text-align:center;"> 2.4895647 </td>
   <td style="text-align:center;"> 0.0438596 </td>
   <td style="text-align:center;"> 0.8598111 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> 26 </td>
   <td style="text-align:center;"> 85.92276 </td>
   <td style="text-align:center;"> 3.7541280 </td>
   <td style="text-align:center;"> 0.3394872 </td>
   <td style="text-align:center;"> 0.0902240 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 57.63293 </td>
   <td style="text-align:center;"> 3.8525660 </td>
   <td style="text-align:center;"> 0.2197802 </td>
   <td style="text-align:center;"> 0.4703526 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 169.70400 </td>
   <td style="text-align:center;"> 4.4133912 </td>
   <td style="text-align:center;"> 0.3956044 </td>
   <td style="text-align:center;"> 0.1821958 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 654.85949 </td>
   <td style="text-align:center;"> 4.3176242 </td>
   <td style="text-align:center;"> -0.3214286 </td>
   <td style="text-align:center;"> 0.4976190 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Juglans regia </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 47.43156 </td>
   <td style="text-align:center;"> 1.0989006 </td>
   <td style="text-align:center;"> -0.0392857 </td>
   <td style="text-align:center;"> 0.8929284 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> 25 </td>
   <td style="text-align:center;"> 17.35603 </td>
   <td style="text-align:center;"> 1.4352301 </td>
   <td style="text-align:center;"> 0.1561538 </td>
   <td style="text-align:center;"> 0.4543071 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 27.89545 </td>
   <td style="text-align:center;"> 3.7114974 </td>
   <td style="text-align:center;"> 0.2120301 </td>
   <td style="text-align:center;"> 0.3678291 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 38.38249 </td>
   <td style="text-align:center;"> 1.4088867 </td>
   <td style="text-align:center;"> -0.2622549 </td>
   <td style="text-align:center;"> 0.3079967 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Mangifera indica </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 17.87202 </td>
   <td style="text-align:center;"> 0.8990063 </td>
   <td style="text-align:center;"> 0.3834586 </td>
   <td style="text-align:center;"> 0.0960014 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 28.80229 </td>
   <td style="text-align:center;"> 0.3127434 </td>
   <td style="text-align:center;"> -0.0392359 </td>
   <td style="text-align:center;"> 0.8771588 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Momordica charantia </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 26.50729 </td>
   <td style="text-align:center;"> 0.9642987 </td>
   <td style="text-align:center;"> 0.0607906 </td>
   <td style="text-align:center;"> 0.8675110 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Nelumbo nucifera </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 97.84551 </td>
   <td style="text-align:center;"> 1.0173700 </td>
   <td style="text-align:center;"> 0.1428571 </td>
   <td style="text-align:center;"> 0.7520337 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 28.16253 </td>
   <td style="text-align:center;"> 1.8005960 </td>
   <td style="text-align:center;"> -0.0489510 </td>
   <td style="text-align:center;"> 0.8863351 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 31.10379 </td>
   <td style="text-align:center;"> 1.2220253 </td>
   <td style="text-align:center;"> 0.4125874 </td>
   <td style="text-align:center;"> 0.1844807 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 56.37928 </td>
   <td style="text-align:center;"> 1.8065166 </td>
   <td style="text-align:center;"> -0.5666667 </td>
   <td style="text-align:center;"> 0.1205743 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 46.80187 </td>
   <td style="text-align:center;"> 3.5497073 </td>
   <td style="text-align:center;"> 0.5363636 </td>
   <td style="text-align:center;"> 0.0936425 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 24.98833 </td>
   <td style="text-align:center;"> 1.1443345 </td>
   <td style="text-align:center;"> -0.3571429 </td>
   <td style="text-align:center;"> 0.4444444 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 28.21185 </td>
   <td style="text-align:center;"> 1.2991090 </td>
   <td style="text-align:center;"> 0.3809524 </td>
   <td style="text-align:center;"> 0.3598710 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Quercus sp </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 59.72765 </td>
   <td style="text-align:center;"> 1.0266447 </td>
   <td style="text-align:center;"> -0.2447552 </td>
   <td style="text-align:center;"> 0.4435043 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Raphanus sativus </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 38.78895 </td>
   <td style="text-align:center;"> 1.2223674 </td>
   <td style="text-align:center;"> 0.7714286 </td>
   <td style="text-align:center;"> 0.1027778 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 16.40985 </td>
   <td style="text-align:center;"> 0.7610006 </td>
   <td style="text-align:center;"> 0.5879121 </td>
   <td style="text-align:center;"> 0.0381070 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 44.58849 </td>
   <td style="text-align:center;"> 2.6615801 </td>
   <td style="text-align:center;"> 0.3666667 </td>
   <td style="text-align:center;"> 0.3362599 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 67.26872 </td>
   <td style="text-align:center;"> 4.3148913 </td>
   <td style="text-align:center;"> 0.5804196 </td>
   <td style="text-align:center;"> 0.0520862 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 60.41812 </td>
   <td style="text-align:center;"> 2.4960328 </td>
   <td style="text-align:center;"> 0.4615385 </td>
   <td style="text-align:center;"> 0.1338363 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 68.36450 </td>
   <td style="text-align:center;"> 2.0348211 </td>
   <td style="text-align:center;"> -0.6484848 </td>
   <td style="text-align:center;"> 0.0490426 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 33.04562 </td>
   <td style="text-align:center;"> 2.0835986 </td>
   <td style="text-align:center;"> 0.1878788 </td>
   <td style="text-align:center;"> 0.6075669 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 667.82536 </td>
   <td style="text-align:center;"> 3.9899830 </td>
   <td style="text-align:center;"> 0.1353383 </td>
   <td style="text-align:center;"> 0.5681081 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum dicoccoides </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 582.40116 </td>
   <td style="text-align:center;"> 4.6545740 </td>
   <td style="text-align:center;"> -0.3608715 </td>
   <td style="text-align:center;"> 0.3056130 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum urartu </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 665.94253 </td>
   <td style="text-align:center;"> 2.3667467 </td>
   <td style="text-align:center;"> 0.0714286 </td>
   <td style="text-align:center;"> 0.9063492 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 43.04148 </td>
   <td style="text-align:center;"> 1.9663878 </td>
   <td style="text-align:center;"> 0.5363636 </td>
   <td style="text-align:center;"> 0.0936425 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 22.43032 </td>
   <td style="text-align:center;"> 1.7051993 </td>
   <td style="text-align:center;"> 0.0175439 </td>
   <td style="text-align:center;"> 0.9454163 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 210.63381 </td>
   <td style="text-align:center;"> 3.5198619 </td>
   <td style="text-align:center;"> 0.0303030 </td>
   <td style="text-align:center;"> 0.9457098 </td>
  </tr>
</tbody>
</table>





```r
S6_peripherybias = species.stats[, c(1:4, 9:10)]
colnames(S6_peripherybias) = c("Species", "Chromosome number", "Chromosome size", 
    "Mean periphery-bias", "Spearman correlation", "p-value")

write.xlsx(x = S6_peripherybias, file = paste(wd, "tables/article_one_supp/tables_supplementary.xls", 
    sep = ""), sheetName = "S6_peripherybias", row.names = FALSE, append = TRUE)
rm(S6_peripherybias)
```



Is the correlation between the periphery-bias ratio and chromosome size significant across species? Chromosomes pooled per species.


```r
nrow(species.stats)
```

```
## [1] 57
```

```r
# Model selection Linear regression
lm_model = lm(peripherybias_ratio ~ phys.map.length, data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = peripherybias_ratio ~ phys.map.length, data = species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.3207 -0.8261 -0.2908  0.6577  2.6611 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     1.746358   0.176310   9.905 7.85e-14 ***
## phys.map.length 0.004417   0.000886   4.985 6.54e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.122 on 55 degrees of freedom
## Multiple R-squared:  0.3112,	Adjusted R-squared:  0.2987 
## F-statistic: 24.85 on 1 and 55 DF,  p-value: 6.536e-06
```

```r
AIC(lm_model)
```

```
## [1] 178.8955
```

```r
BIC(lm_model)
```

```
## [1] 185.0246
```

```r
# Quadratic regression
lm_model = lm(peripherybias_ratio ~ phys.map.length + I(phys.map.length^2), data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = peripherybias_ratio ~ phys.map.length + I(phys.map.length^2), 
##     data = species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.8388 -0.6438 -0.2247  0.5360  2.0807 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           1.158e+00  2.226e-01   5.203 3.11e-06 ***
## phys.map.length       1.751e-02  3.573e-03   4.901 9.09e-06 ***
## I(phys.map.length^2) -2.052e-05  5.458e-06  -3.760  0.00042 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.009 on 54 degrees of freedom
## Multiple R-squared:  0.4541,	Adjusted R-squared:  0.4339 
## F-statistic: 22.46 on 2 and 54 DF,  p-value: 7.985e-08
```

```r
AIC(lm_model)
```

```
## [1] 167.6433
```

```r
BIC(lm_model)
```

```
## [1] 175.8155
```

Selected model for simplicity: LM


```r
# The mean species correlation estimated by 1,000 bootstraps
mean(boot)
```

```
## [1] 0.1201095
```

```r
quantile(boot, 0.025)
```

```
##       2.5% 
## 0.01758947
```

```r
quantile(boot, 0.975)
```

```
##     97.5% 
## 0.2195359
```

```r
# Species correlations
summary(species.stats$correlation_peripherybias)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -1.00000 -0.07143  0.13534  0.12078  0.39560  0.80952
```

```r
cor.test(species.stats$peripherybias_ratio, species.stats$phys.map.length, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  species.stats$peripherybias_ratio and species.stats$phys.map.length
## S = 12226, p-value = 1.175e-06
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.6037724
```

```r
# Linear regression
lm_model = lm(peripherybias_ratio ~ log10(phys.map.length), data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = peripherybias_ratio ~ log10(phys.map.length), data = species.stats)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -2.24893 -0.47523 -0.09788  0.53007  2.06181 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)             -1.2315     0.5314  -2.317   0.0242 *  
## log10(phys.map.length)   1.9931     0.2972   6.705 1.15e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.003 on 55 degrees of freedom
## Multiple R-squared:  0.4498,	Adjusted R-squared:  0.4398 
## F-statistic: 44.96 on 1 and 55 DF,  p-value: 1.145e-08
```



```
## [1] 2.190423
```

```
##     2.5% 
## 2.062457
```

```
##    97.5% 
## 2.322344
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-82-1.png" alt="(ref:peripherybias-species)"  />
<p class="caption">(\#fig:unnamed-chunk-82-1)(ref:peripherybias-species)</p>
</div><div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-82-2.png" alt="(ref:peripherybias-species)"  />
<p class="caption">(\#fig:unnamed-chunk-82-2)(ref:peripherybias-species)</p>
</div>

(ref:peripherybias-species) (ref:peripherybias-species) The periphery-bias ratio depends on the mean chromosome size across species (n = 38). The linear regression line and its parametric 95% confidence interval were estimated in ggplot2. A periphery-bias ratio above one indicates that recombination rate in the tips of the chromosome are higher than the mean chromosome recombination rate. The inset presents the distribution of periphery-bias ratios (n = 665 chromosomes). The mean periphery-bias ratio and its 95% confidence interval (black solid and dashed lines) were estimated by 1,000 bootstraps. The theoretical value for equal recombination in the tips than the rest of the chromosome (periphery-bias ratio of one) is the red vertical line.







```r
lm.mod = lm(peripherybias_ratio ~ linkage.map.length, data = tmp_species.stats)
summary(lm.mod)
```

```
## 
## Call:
## lm(formula = peripherybias_ratio ~ linkage.map.length, data = tmp_species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.9444 -1.0036 -0.3173  0.9565  3.3716 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)   
## (Intercept)        1.780529   0.543201   3.278  0.00182 **
## linkage.map.length 0.003550   0.004158   0.854  0.39688   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.344 on 55 degrees of freedom
## Multiple R-squared:  0.01308,	Adjusted R-squared:  -0.004861 
## F-statistic: 0.7291 on 1 and 55 DF,  p-value: 0.3969
```


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-85-1.png" alt="(ref:peripherybias-linkagemaplength)"  />
<p class="caption">(\#fig:unnamed-chunk-85)(ref:peripherybias-linkagemaplength)</p>
</div>
(ref:peripherybias-linkagemaplength) Linkage map length (i.e. crossover frequency) does not correlate with the periphery-bias ratio.







Is the relationship between chromosome size and the periphery-bias ratio universal within species?


```r
qqPlot(sqrt(chromosome.stats$peripherybias_ratio), ylab = "Periphery-bias ratio")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-86-1.png)<!-- -->

```
## [1] 103 310
```

```r
cor.test(log10(chromosome.stats$phys.map.length), sqrt(chromosome.stats$peripherybias_ratio), 
    method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  log10(chromosome.stats$phys.map.length) and sqrt(chromosome.stats$peripherybias_ratio)
## S = 23534674, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.4928681
```

```r
# LINEAR MODEL
lm_model = lm(sqrt(peripherybias_ratio) ~ log10(phys.map.length), data = chromosome.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = sqrt(peripherybias_ratio) ~ log10(phys.map.length), 
##     data = chromosome.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.3270 -0.3129 -0.0238  0.3091  1.3359 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)             0.27550    0.07516   3.666 0.000267 ***
## log10(phys.map.length)  0.64989    0.04333  15.000  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.4927 on 651 degrees of freedom
##   (12 observations deleted due to missingness)
## Multiple R-squared:  0.2568,	Adjusted R-squared:  0.2557 
## F-statistic:   225 on 1 and 651 DF,  p-value: < 2.2e-16
```

```r
AIC(lm_model)
```

```
## [1] 932.7827
```

```r
# LINEAR MIXED MODEL
lmer_model = lmer(sqrt(peripherybias_ratio) ~ log10(phys.map.length) + (1 | species), 
    data = chromosome.stats)
summary(lmer_model)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: sqrt(peripherybias_ratio) ~ log10(phys.map.length) + (1 | species)
##    Data: chromosome.stats
## 
## REML criterion at convergence: 726.4
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.8023 -0.5813  0.0342  0.5713  3.0047 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  species  (Intercept) 0.09366  0.3060  
##  Residual             0.14790  0.3846  
## Number of obs: 653, groups:  species, 57
## 
## Fixed effects:
##                         Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)              0.25031    0.14938 102.50053   1.676   0.0968 .  
## log10(phys.map.length)   0.65566    0.08353 112.44430   7.849 2.67e-12 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## lg10(phy..) -0.956
```

```r
AIC(lmer_model)
```

```
## [1] 734.4097
```


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-87-1.png" alt="(ref:chromosomesize-peripherybiasratio)"  />
<p class="caption">(\#fig:unnamed-chunk-87)(ref:chromosomesize-peripherybiasratio)</p>
</div>

(ref:chromosomesize-peripherybiasratio) The periphery-bias ratio increases with chromosome size across species, but the relationship does not stand within species.


### Sensitivity of the metric (periphery-bias ratio) to the sampling scale.


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-88-1.png" alt="(ref:peripherybias-sensitivity)"  />
<p class="caption">(\#fig:unnamed-chunk-88)(ref:peripherybias-sensitivity)</p>
</div>

(ref:peripherybias-sensitivity) Sensitivity of the periphery-bias ratio metric to the sampling scale (i.e. the number of bins sampled at the tip of the chromosome). The periphery-bias ratio seems steady with an increasing number of bins sampled, though the ratio is slightly decreasing when more bins are sampled.


#### Figure S12





# Telomeres and centromeres both influence the distribution of crossovers


## Telomere effect

Recombination rates increase with the distance to the nearest telomere.

Correlations between recombination rates and distance to the nearest telomere.

Test which species have a significant correlation with the distance to the telomere (supplementary table).





```r
S7_correlationchromosomes = cor_dist2telomere[, 1:4]
colnames(S7_correlationchromosomes) = c("Species", "Chromosome", "Spearman correlation", 
    "p-value")

write.xlsx(x = S7_correlationchromosomes, file = paste(wd, "tables/article_one_supp/tables_supplementary.xls", 
    sep = ""), sheetName = "S7_correlationchromosomes", row.names = FALSE, append = TRUE)
rm(S7_correlationchromosomes)
```



Correlations at species and chromosome level




```r
# Sample size
nrow(species.stats)
```

```
## [1] 57
```



```r
kable(species.stats[, c(1:3, 11)], caption = "Species correlation between the distance to the nearest telomere and the recombination rate.", 
    labels = c("Species", "Chromosome number", "Chromosome size", "Spearman correlation"), 
    align = "c")
```

<table>
<caption>(\#tab:unnamed-chunk-94)Species correlation between the distance to the nearest telomere and the recombination rate.</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> species </th>
   <th style="text-align:center;"> chromosome_number </th>
   <th style="text-align:center;"> phys.map.length </th>
   <th style="text-align:center;"> correlation_dist2telomere </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 517.44324 </td>
   <td style="text-align:center;"> -0.9447383 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 23.82927 </td>
   <td style="text-align:center;"> -0.1925520 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 100.95097 </td>
   <td style="text-align:center;"> -0.5861152 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 124.68497 </td>
   <td style="text-align:center;"> -0.7865345 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Boechera stricta </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 26.52555 </td>
   <td style="text-align:center;"> -0.4861575 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 45.62172 </td>
   <td style="text-align:center;"> -0.7064544 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 44.71564 </td>
   <td style="text-align:center;"> -0.3911250 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 26.14066 </td>
   <td style="text-align:center;"> -0.5416268 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 30.63062 </td>
   <td style="text-align:center;"> -0.3028815 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 222.10960 </td>
   <td style="text-align:center;"> -0.0965251 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 15.83906 </td>
   <td style="text-align:center;"> 0.5766299 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsicum annuum </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 227.12175 </td>
   <td style="text-align:center;"> -0.8177649 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cenchrus americanus </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 223.50536 </td>
   <td style="text-align:center;"> -0.3117457 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrullus lanatus </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 33.59728 </td>
   <td style="text-align:center;"> -0.4349560 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 26.02017 </td>
   <td style="text-align:center;"> -0.6789178 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Coffea canephora </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 53.00321 </td>
   <td style="text-align:center;"> -0.7962224 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 30.99121 </td>
   <td style="text-align:center;"> -0.7843992 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 27.40843 </td>
   <td style="text-align:center;"> -0.0274564 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 10.34886 </td>
   <td style="text-align:center;"> -0.2031031 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 10.78700 </td>
   <td style="text-align:center;"> -0.1623203 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 23.84612 </td>
   <td style="text-align:center;"> -0.7774816 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Draba nivalis </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 34.11541 </td>
   <td style="text-align:center;"> -0.7419194 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 46.94000 </td>
   <td style="text-align:center;"> 0.2496445 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 55.68998 </td>
   <td style="text-align:center;"> -0.3969593 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 48.12719 </td>
   <td style="text-align:center;"> -0.6771574 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> 26 </td>
   <td style="text-align:center;"> 85.92276 </td>
   <td style="text-align:center;"> -0.8310385 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 57.63293 </td>
   <td style="text-align:center;"> -0.7418663 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 169.70400 </td>
   <td style="text-align:center;"> -0.6588514 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 654.85949 </td>
   <td style="text-align:center;"> -0.8248608 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Juglans regia </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 47.43156 </td>
   <td style="text-align:center;"> -0.4387363 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> 25 </td>
   <td style="text-align:center;"> 17.35603 </td>
   <td style="text-align:center;"> -0.5677997 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 27.89545 </td>
   <td style="text-align:center;"> -0.7649821 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 38.38249 </td>
   <td style="text-align:center;"> -0.6837059 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Mangifera indica </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 17.87202 </td>
   <td style="text-align:center;"> -0.2076802 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 28.80229 </td>
   <td style="text-align:center;"> -0.0641793 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Momordica charantia </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 26.50729 </td>
   <td style="text-align:center;"> -0.4024631 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Nelumbo nucifera </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 97.84551 </td>
   <td style="text-align:center;"> 0.1356657 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 28.16253 </td>
   <td style="text-align:center;"> -0.5404312 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 31.10379 </td>
   <td style="text-align:center;"> -0.5798772 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 56.37928 </td>
   <td style="text-align:center;"> -0.8367877 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 46.80187 </td>
   <td style="text-align:center;"> -0.6929217 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 24.98833 </td>
   <td style="text-align:center;"> -0.1507408 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 28.21185 </td>
   <td style="text-align:center;"> -0.6062939 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Quercus sp </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 59.72765 </td>
   <td style="text-align:center;"> -0.0568060 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Raphanus sativus </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 38.78895 </td>
   <td style="text-align:center;"> -0.1157651 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 16.40985 </td>
   <td style="text-align:center;"> -0.0026683 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 44.58849 </td>
   <td style="text-align:center;"> -0.7556303 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 67.26872 </td>
   <td style="text-align:center;"> -0.7020987 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 60.41812 </td>
   <td style="text-align:center;"> -0.6201034 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 68.36450 </td>
   <td style="text-align:center;"> -0.7765123 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 33.04562 </td>
   <td style="text-align:center;"> -0.7852950 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 667.82536 </td>
   <td style="text-align:center;"> -0.8434463 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum dicoccoides </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 582.40116 </td>
   <td style="text-align:center;"> -0.8247954 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum urartu </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 665.94253 </td>
   <td style="text-align:center;"> -0.8544557 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 43.04148 </td>
   <td style="text-align:center;"> -0.6674743 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 22.43032 </td>
   <td style="text-align:center;"> -0.4937258 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 210.63381 </td>
   <td style="text-align:center;"> -0.8078647 </td>
  </tr>
</tbody>
</table>



```r
S8_correlationspecies = species.stats[, c(1:3, 11)]
colnames(S8_correlationspecies) = c("Species", "Chromosome number", "Chromosome size", 
    "Spearman correlation")

write.xlsx(x = S8_correlationspecies, file = paste(wd, "tables/article_one_supp/tables_supplementary.xls", 
    sep = ""), sheetName = "S8_correlationspecies", row.names = FALSE, append = TRUE)
rm(S8_correlationspecies)
```



Selected model for simplicity: LM


```r
# The mean species correlation estimated by 1,000 bootstraps
mean(boot)
```

```
## [1] -0.5044995
```

```r
quantile(boot, 0.025)
```

```
##       2.5% 
## -0.5868395
```

```r
quantile(boot, 0.975)
```

```
##      97.5% 
## -0.4142717
```

```r
# Species correlations
summary(species.stats$correlation_dist2telomere)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -0.9447 -0.7765 -0.6063 -0.5050 -0.3029  0.5766
```

```r
cor.test(species.stats$correlation_dist2telomere, species.stats$phys.map.length, 
    method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  species.stats$correlation_dist2telomere and species.stats$phys.map.length
## S = 46712, p-value = 5.565e-05
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.5138709
```

```r
# Linear regression
lm_model = lm(correlation_dist2telomere ~ phys.map.length, data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = correlation_dist2telomere ~ phys.map.length, data = species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.3673 -0.2302 -0.0907  0.2341  1.0177 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     -0.4299957  0.0485244  -8.861 3.52e-12 ***
## phys.map.length -0.0007009  0.0002438  -2.874  0.00574 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3089 on 55 degrees of freedom
## Multiple R-squared:  0.1306,	Adjusted R-squared:  0.1148 
## F-statistic: 8.263 on 1 and 55 DF,  p-value: 0.005744
```



<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-97-1.png" alt="(ref:dist2telomere-species)"  />
<p class="caption">(\#fig:unnamed-chunk-97)(ref:dist2telomere-species)</p>
</div>

(ref:dist2telomere-species) (ref:dist2telomere-species) The negative correlation between recombination rates and the distance to the nearest telomere is stronger for species with a larger chromosome size (n = 57). The linear regression line and its parametric 95% confidence interval were estimated in ggplot2. The inset presents the distribution of Spearman correlation coefficients for chromosomes (n = 665 chromosomes). The mean correlation and its 95% confidence interval (black solid and dashed lines) were estimated by 1,000 bootstraps. The red vertical line is for a null correlation.

#### Figure S5




Interestingly, at a species level, we can identify two groups of species. Although most species correlations are globally significant (i.e. a majority of significant chromosomes), there were a group of high negative correlations (lower than -0.5) and another one of lower/null correlations (between -0.5 and 0.1).

<table>
<caption>(\#tab:SpeciesCorrelationPercentage of significant correlationsChromosome length)Species correlation between chromosome size and distance to the nearest telomere</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> species </th>
   <th style="text-align:center;"> correlation </th>
   <th style="text-align:center;"> signif </th>
   <th style="text-align:center;"> phys.map.length </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> -0.94 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 517.44 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> -0.19 </td>
   <td style="text-align:center;"> 0.40 </td>
   <td style="text-align:center;"> 23.83 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> -0.59 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 100.95 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> -0.79 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 124.68 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Boechera stricta </td>
   <td style="text-align:center;"> -0.49 </td>
   <td style="text-align:center;"> 0.71 </td>
   <td style="text-align:center;"> 26.53 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> -0.71 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 45.62 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> -0.39 </td>
   <td style="text-align:center;"> 0.83 </td>
   <td style="text-align:center;"> 44.72 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> -0.54 </td>
   <td style="text-align:center;"> 0.86 </td>
   <td style="text-align:center;"> 26.14 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> -0.30 </td>
   <td style="text-align:center;"> 0.87 </td>
   <td style="text-align:center;"> 30.63 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> -0.10 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 222.11 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> 0.58 </td>
   <td style="text-align:center;"> 0.86 </td>
   <td style="text-align:center;"> 15.84 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsicum annuum </td>
   <td style="text-align:center;"> -0.82 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 227.12 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cenchrus americanus </td>
   <td style="text-align:center;"> -0.31 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 223.51 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrullus lanatus </td>
   <td style="text-align:center;"> -0.43 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 33.60 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> -0.68 </td>
   <td style="text-align:center;"> 0.88 </td>
   <td style="text-align:center;"> 26.02 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Coffea canephora </td>
   <td style="text-align:center;"> -0.80 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 53.00 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> -0.78 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 30.99 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> -0.03 </td>
   <td style="text-align:center;"> 0.86 </td>
   <td style="text-align:center;"> 27.41 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> -0.20 </td>
   <td style="text-align:center;"> 0.58 </td>
   <td style="text-align:center;"> 10.35 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> -0.16 </td>
   <td style="text-align:center;"> 0.71 </td>
   <td style="text-align:center;"> 10.79 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> -0.78 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 23.85 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Draba nivalis </td>
   <td style="text-align:center;"> -0.74 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 34.12 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> 0.25 </td>
   <td style="text-align:center;"> 0.82 </td>
   <td style="text-align:center;"> 46.94 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> -0.40 </td>
   <td style="text-align:center;"> 0.82 </td>
   <td style="text-align:center;"> 55.69 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> -0.68 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 48.13 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> -0.83 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 85.92 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> -0.74 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 57.63 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> -0.66 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 169.70 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> -0.82 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 654.86 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Juglans regia </td>
   <td style="text-align:center;"> -0.44 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 47.43 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> -0.57 </td>
   <td style="text-align:center;"> 0.92 </td>
   <td style="text-align:center;"> 17.36 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> -0.76 </td>
   <td style="text-align:center;"> 0.95 </td>
   <td style="text-align:center;"> 27.90 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> -0.68 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 38.38 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Mangifera indica </td>
   <td style="text-align:center;"> -0.21 </td>
   <td style="text-align:center;"> 0.70 </td>
   <td style="text-align:center;"> 17.87 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> -0.06 </td>
   <td style="text-align:center;"> 0.83 </td>
   <td style="text-align:center;"> 28.80 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Momordica charantia </td>
   <td style="text-align:center;"> -0.40 </td>
   <td style="text-align:center;"> 0.82 </td>
   <td style="text-align:center;"> 26.51 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Nelumbo nucifera </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 97.85 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> -0.54 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 28.16 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> -0.58 </td>
   <td style="text-align:center;"> 0.83 </td>
   <td style="text-align:center;"> 31.10 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> -0.84 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 56.38 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> -0.69 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 46.80 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> -0.15 </td>
   <td style="text-align:center;"> 0.43 </td>
   <td style="text-align:center;"> 24.99 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> -0.61 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 28.21 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Quercus sp </td>
   <td style="text-align:center;"> -0.06 </td>
   <td style="text-align:center;"> 0.83 </td>
   <td style="text-align:center;"> 59.73 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Raphanus sativus </td>
   <td style="text-align:center;"> -0.12 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 38.79 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> 0.00 </td>
   <td style="text-align:center;"> 0.85 </td>
   <td style="text-align:center;"> 16.41 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> -0.76 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 44.59 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> -0.70 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 67.27 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> -0.62 </td>
   <td style="text-align:center;"> 0.92 </td>
   <td style="text-align:center;"> 60.42 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> -0.78 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 68.36 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> -0.79 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 33.05 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> -0.84 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 667.83 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum dicoccoides </td>
   <td style="text-align:center;"> -0.82 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 582.40 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum urartu </td>
   <td style="text-align:center;"> -0.85 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 665.94 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> -0.67 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 43.04 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> -0.49 </td>
   <td style="text-align:center;"> 0.84 </td>
   <td style="text-align:center;"> 22.43 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> -0.81 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 210.63 </td>
  </tr>
</tbody>
</table>


Finally, to identify an inter-specific recurrent pattern, chromosomes were pooled per species. The relationship between recombination rates and the distance toward telomeres was assessed by quadratic regression with a Liner Mixed Model specifying a random species effect. Coefficients of the model and their confidence intervals at 95% were estimated with parametric bootstrap.





<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-101-1.png" alt="(ref:DistancesRelative-pooledChromosomes)"  />
<p class="caption">(\#fig:unnamed-chunk-101)(ref:DistancesRelative-pooledChromosomes)</p>
</div>

(ref:DistancesRelative-pooledChromosomes) Standardized recombination rate (cM/Mb) as a function of the relative distance from the telomere along the chromosome (physical distances expressed in 20 bins). Chromosomes were split in halves, distance 0.5 being the center of the chromosome. Then, chromosomes were pooled per species. Each color is a species. Linear quadratic regression estimated with a Linear Model is the black line (95% confidence interval in dark grey).


To uncover specific effects that could ultimately differ from the quadratic regression, I plotted loess regression lines over the pooled data.

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-102-1.png" alt="(ref:DistancesRelative-pooledChromosomes-LOESS)"  />
<p class="caption">(\#fig:unnamed-chunk-102)(ref:DistancesRelative-pooledChromosomes-LOESS)</p>
</div>

(ref:DistancesRelative-pooledChromosomes-LOESS) Standardized recombination rate (cM/Mb) as a function of the relative distance from the telomere along the chromosome (physical distances expressed in 20 bins). Chromosomes were split in halves, distance 0.5 being the center of the chromosome. Then, chromosomes were pooled per species. Each color is a species. A loess regression was estimated for each species.


<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-103-1.png" style="display: block; margin: auto;" /><img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-103-2.png" style="display: block; margin: auto;" /><img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-103-3.png" style="display: block; margin: auto;" /><img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-103-4.png" style="display: block; margin: auto;" /><img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-103-5.png" style="display: block; margin: auto;" />

#### Figure S6





Classify species in three categories: telomere pattern, sub-telomere pattern and exceptions.




Exception patterns (6 species) tend to be higher at the center than on the rest of the chromosome. Interestingly, *Nelumbo nucifera* and *Camellia sinensis* are among these outlying patterns.


```r
kable(telomere_pattern)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> species </th>
   <th style="text-align:left;"> pattern </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Aegilops speltoides </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arabidopsis thaliana </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arachis duranensis </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Arachis hypogaea </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Boechera stricta </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Brachypodium distachyon </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Brassica napus </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Brassica rapa </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Camelina sativa </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Camellia sinensis </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Capsella rubella </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Capsicum annuum </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cenchrus americanus </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Citrullus lanatus </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Citrus sinensis </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Coffea canephora </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cucumis melo </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cucumis sativus </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cucurbita maxima </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cucurbita pepo </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Dioscorea alata </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Draba nivalis </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Elaeis guineensis </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Eucalyptus grandis </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Glycine max </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gossypium hirsutum </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gossypium raimondii </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Helianthus annuus </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Hordeum vulgare </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Juglans regia </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lupinus albus </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lupinus angustifolius </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Malus domestica </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mangifera indica </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Manihot esculenta </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Momordica charantia </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nelumbo nucifera </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oryza nivara </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Oryza sativa </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Panicum hallii </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Phaseolus vulgaris </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Prunus mume </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Prunus persica </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Quercus sp </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Raphanus sativus </td>
   <td style="text-align:left;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sesamum indicum </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Setaria italica </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Solanum lycopersicum </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Solanum tuberosum </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sorghum bicolor </td>
   <td style="text-align:left;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Theobroma cacao </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Triticum aestivum </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Triticum dicoccoides </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Triticum urartu </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vigna unguiculata </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Vitis vinifera </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Zea mays </td>
   <td style="text-align:left;"> telomere </td>
  </tr>
</tbody>
</table>



```r
# The distribution of patterns among species
table(telomere_pattern$pattern)
```

```
## 
## exception      peak  telomere 
##         7        16        34
```


```r
chromosome.stats.telomere = merge(chromosome.stats, telomere_pattern)
boxplot(phys.map.length ~ pattern, data = chromosome.stats.telomere)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-108-1.png)<!-- -->

```r
summary(chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == 
    "telomere"], na.rm = TRUE)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   10.68   28.64   48.58  128.38  115.51  830.83
```

```r
summary(chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == 
    "peak"], na.rm = TRUE)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    5.05   17.22   27.22   29.88   37.14   80.88
```

```r
summary(chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == 
    "exception"], na.rm = TRUE)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   13.32   31.61   51.70   68.25   72.22  336.97
```

```r
wilcox.test(chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == 
    "telomere"], chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == 
    "peak"])
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == "telomere"] and chromosome.stats.telomere$phys.map.length[chromosome.stats.telomere$pattern == "peak"]
## W = 61080, p-value < 2.2e-16
## alternative hypothesis: true location shift is not equal to 0
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-109-1.png" alt="(ref:DistancesRelative-pooledChromosomes-patterns)"  />
<p class="caption">(\#fig:unnamed-chunk-109)(ref:DistancesRelative-pooledChromosomes-patterns)</p>
</div>


(ref:DistancesRelative-pooledChromosomes-patterns) (ref:DistancesRelative-pooledChromosomes-patterns) Recombination rates are higher in distal regions and lower near the center of the chromosome. Standardized recombination rates (cM/Mb) as a function of the relative distance from the telomere along the chromosome (physical distances expressed in 20 bins). Two patterns were identified and species were pooled accordingly. Recombination rates decrease immediately from the tip of the chromosome for 35 species (dark grey line and ribbon) or recombination rate is reduced in telomeric regions and the peak of recombination is in a sub-telomeric region for 16 species (light grey line and ribbon). The solid line represents the mean recombination rate estimated in a bin and upper and lower boundaries of the ribbon represent the maximum and minimum values for a particular pattern. Besides, patterns that were not classified (6 species) were represented by loess regression in grey dashed lines. For estimating recombination rates in bins of relative distance, chromosomes were split in halves, a distance of 0.5 being the center of the chromosome. Then, chromosomes were pooled per species.

#### Figure 6




### Figure 7






Redraw the periphery-bias ratio with colors according to the classification of patterns.


```r
# The mean species correlation estimated by 1,000 bootstraps
mean(boot)
```

```
## [1] -0.5044995
```

```r
quantile(boot, 0.025)
```

```
##       2.5% 
## -0.5868395
```

```r
quantile(boot, 0.975)
```

```
##      97.5% 
## -0.4142717
```

```r
# Species correlations
summary(species.stats$correlation_peripherybias)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -1.00000 -0.07143  0.13534  0.12078  0.39560  0.80952
```

```r
cor.test(species.stats$peripherybias_ratio, species.stats$phys.map.length, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  species.stats$peripherybias_ratio and species.stats$phys.map.length
## S = 12226, p-value = 1.175e-06
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.6037724
```

```r
# Linear regression
lm_model = lm(peripherybias_ratio ~ phys.map.length, data = species.stats)
summary(lm_model)
```

```
## 
## Call:
## lm(formula = peripherybias_ratio ~ phys.map.length, data = species.stats)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -2.3207 -0.8261 -0.2908  0.6577  2.6611 
## 
## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     1.746358   0.176310   9.905 7.85e-14 ***
## phys.map.length 0.004417   0.000886   4.985 6.54e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.122 on 55 degrees of freedom
## Multiple R-squared:  0.3112,	Adjusted R-squared:  0.2987 
## F-statistic: 24.85 on 1 and 55 DF,  p-value: 6.536e-06
```


```
## [1] 2.192275
```

```
##     2.5% 
## 2.066465
```

```
##    97.5% 
## 2.323914
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-113-1.png" alt="(ref:peripherybias-species)"  />
<p class="caption">(\#fig:unnamed-chunk-113-1)(ref:peripherybias-species)</p>
</div><div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-113-2.png" alt="(ref:peripherybias-species)"  />
<p class="caption">(\#fig:unnamed-chunk-113-2)(ref:peripherybias-species)</p>
</div>

(ref:peripherybias-species) (ref:peripherybias-species) The periphery-bias ratio depends on the mean chromosome size across species (n = 38). The linear regression line and its parametric 95% confidence interval were estimated in ggplot2. A periphery-bias ratio above one indicates that recombination rate in the tips of the chromosome are higher than the mean chromosome recombination rate. The inset presents the distribution of periphery-bias ratios (n = 665 chromosomes). The mean periphery-bias ratio and its 95% confidence interval (black solid and dashed lines) were estimated by 1,000 bootstraps. The theoretical value for equal recombination in the tips than the rest of the chromosome (periphery-bias ratio of one) is the red vertical line.

#### Figure 5






## Centromere effect

### How many chromosomes have lower recombination in the centromere?

For each chromosome, estimate the recombination rate within the centromere (recombination rate of the 100 kb window overlapping the centromere position).





```r
# How many chromosomes with a centromere position?
sum(!is.na(chromosome.stats$rec.in.centromere))
```

```
## [1] 425
```

```r
# How many chromosomes with no recombination in centromere
sum(chromosome.stats$rec.in.centromere == 0, na.rm = TRUE)
```

```
## [1] 81
```

```r
sum(chromosome.stats$rec.in.centromere == 0, na.rm = TRUE)/sum(!is.na(chromosome.stats$rec.in.centromere))
```

```
## [1] 0.1905882
```

```r
# How many chromosomes with less recombination in centromere than the chromosome
# averaged recombination rate?
sum(chromosome.stats$rec.in.centromere < chromosome.stats$mean.recrate, na.rm = TRUE)
```

```
## [1] 388
```

```r
sum(chromosome.stats$rec.in.centromere < chromosome.stats$mean.recrate, na.rm = TRUE)/sum(!is.na(chromosome.stats$rec.in.centromere))
```

```
## [1] 0.9129412
```



Resampling test for each chromosome. Is the centromeric rec. rate significantly lower than the chromosome averaged rec. rate?





```r
# How many chromosomes with SIGNIFICANTLY less recombination in centromere than
# the chromosome averaged recombination rate?

# Resampling test. 1000 bootstraps.
sum(chromosome.stats$rec.in.centromere < chromosome.stats$lower.average.recrate, 
    na.rm = TRUE)
```

```
## [1] 385
```

```r
sum(chromosome.stats$rec.in.centromere < chromosome.stats$lower.average.recrate, 
    na.rm = TRUE)/sum(!is.na(chromosome.stats$rec.in.centromere))
```

```
## [1] 0.9058824
```

```r
# And significantly higher?
sum(chromosome.stats$rec.in.centromere > chromosome.stats$upper.average.recrate, 
    na.rm = TRUE)
```

```
## [1] 30
```

```r
sum(chromosome.stats$rec.in.centromere > chromosome.stats$upper.average.recrate, 
    na.rm = TRUE)/sum(!is.na(chromosome.stats$rec.in.centromere))
```

```
## [1] 0.07058824
```



```r
par(mfrow = c(1, 2))
hist(chromosome.stats$rec.in.centromere, main = "", xlab = "Recombination in centromere", 
    breaks = 50)
hist(chromosome.stats$boot.average.recrate, main = "", xlab = "Averaged recombination", 
    breaks = 50)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-119-1.png)<!-- -->

```r
par(mfrow = c(1, 1))
```

### Formal test of a centromere effect

Proportion of chromosomes that are metacentric or sub-metacentric in the dataset (long arm/short arm ratio < 3).


```r
ratio = rep(NA, nrow(chromosome.stats))
for (i in 1:nrow(chromosome.stats)) {
    ratio[i] = max(c(chromosome.stats$centromeric_index[i], (1 - chromosome.stats$centromeric_index[i])))/min(c(chromosome.stats$centromeric_index[i], 
        (1 - chromosome.stats$centromeric_index[i])))
}
hist(ratio, main = "", breaks = 30)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-120-1.png)<!-- -->

```r
sum(ratio < 3, na.rm = TRUE)/sum(!is.na(ratio))
```

```
## [1] 0.92
```

With d(x) the genetic distance at the relative physical position x (x betwwen 0 and 1) and c the centromeric index.

(1) Haenel et al. (2018) telomere model assumes that d(1/2) = d(1) - d(1/2), hence d(1/2)/d(1) = 1/2

(2) If a centromere effect and 1 co/chromosome, d(c) 	= c x 50 + a x c
 and d(1) – d(c) = (1-c) x 50 + a x (1-c). So, d(c) / d(1) = c
 
(3) If a centromere effect and 1 co/chromosome arm, d(c) 	= 50 + a x c
 and d(1) – d(c) = 50 + a x (1-c). So, (d(c) -50) / (d(1) – 100) = c

For chromosomes with a centromeric index, we used the Marey map to compute the predictions of the model.



Verify the predictions.


```r
# How many predictions are negative for one CO per chromosome arm (i.e. shorter
# arm < 50cM)?  Count:
sum((predictions_centromereModels$pred_centromere_arm < 0), na.rm = TRUE)
```

```
## [1] 94
```

```r
# Proportion:
sum((predictions_centromereModels$pred_centromere_arm < 0)/nrow(predictions_centromereModels), 
    na.rm = TRUE)
```

```
## [1] 0.2211765
```




```r
par(mfrow = c(1, 2))
hist(predictions_centromereModels$pred_telomere, xlab = "Predicted ratio", breaks = 30, 
    main = "")
# Distribution centered on 0.5, but too much variance than expected under the
# strict telomere model of Haenel et al. (2018)
plot(rep(0.5, nrow(predictions_centromereModels)), predictions_centromereModels$pred_telomere, 
    xlab = "Predicted", ylab = "Observed", main = "")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-123-1.png)<!-- -->

```r
# The prediction of a ratio constant at 0.5 is rejected

hist(predictions_centromereModels$pred_centromere_chromosome, xlab = "Predicted ratio", 
    breaks = 30, main = "")
plot(predictions_centromereModels$centromeric_index, predictions_centromereModels$pred_centromere_chromosome, 
    xlab = "Predicted", ylab = "Observed", main = "")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-123-2.png)<!-- -->

```r
hist(predictions_centromereModels$pred_centromere_arm, xlab = "Predicted ratio", 
    breaks = 30, main = "")
plot(predictions_centromereModels$centromeric_index, predictions_centromereModels$pred_centromere_arm, 
    xlab = "Predicted", ylab = "Observed", main = "")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-123-3.png)<!-- -->

```r
# Manually remove predictions outside (0,1)
sum(predictions_centromereModels$pred_centromere_arm > 1 | predictions_centromereModels$pred_centromere_arm < 
    0, na.rm = TRUE)/sum(!is.na(predictions_centromereModels$pred_centromere_arm), 
    na.rm = TRUE)
```

```
## [1] 0.4882075
```

```r
hist(predictions_centromereModels$pred_centromere_arm[which(predictions_centromereModels$pred_centromere_arm < 
    1 & predictions_centromereModels$pred_centromere_arm > 0)], xlab = "Predicted ratio", 
    breaks = 30, main = "")
plot(predictions_centromereModels$centromeric_index[which(predictions_centromereModels$pred_centromere_arm < 
    1 & predictions_centromereModels$pred_centromere_arm > 0)], predictions_centromereModels$pred_centromere_arm[which(predictions_centromereModels$pred_centromere_arm < 
    1 & predictions_centromereModels$pred_centromere_arm > 0)], xlab = "Predicted", 
    ylab = "Observed", main = "")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-123-4.png)<!-- -->

```r
par(mfrow = c(1, 1))
```

Estimate the correlation expected-observed for the best model: centromere with one crossover per chromosome.


```r
# Sample size
nrow(predictions_centromereModels)
```

```
## [1] 425
```


```r
cor.test(predictions_centromereModels$centromeric_index, predictions_centromereModels$pred_centromere_chromosome, 
    method = "spearman")
```

```
## Warning in cor.test.default(predictions_centromereModels$centromeric_index, :
## Cannot compute exact p-value with ties
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  predictions_centromereModels$centromeric_index and predictions_centromereModels$pred_centromere_chromosome
## S = 3626920, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.7165184
```

Make the figure.


```
## [1] 0.4539616
```

```
##      2.5% 
## 0.4338565
```

```
##     97.5% 
## 0.4730382
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-126-1.png" alt="(ref:centromere-predictions)"  />
<p class="caption">(\#fig:unnamed-chunk-126)(ref:centromere-predictions)</p>
</div>

<!-- (ref:centromere-predictions) (ref:centromere-predictions) (A) Predictions telomere model. (B) Predictions centromere model with one mandatory crossover. (C) Correlation observed-predicted under the centromere model with one mandatory crossover. -->
(ref:centromere-predictions) (ref:centromere-predictions) Correlation observed-predicted under the centromere model with one mandatory crossover.


### Model selection

Model selection based on R2 across species.


```r
centromere_modselection = data.frame(Model = c("Telomere", "Centromere (arm)", "Centromere (chromosome)", 
    "Telomere (subset)", "Centromere (arm) (subset)", "Centromere (chromosome) (subset)", 
    "Telomere (trueCI)", "Centromere (arm) (trueCI)", "Centromere (chromosome) (trueCI)"), 
    Expected = rep(c("d(1/2) / d(1) = 0.5", "(d(c) -50) / (d(1) – 100) = c", "d(c) / d(1) = c"), 
        3), AdjustedRsquared = NA, Pvalue = NA, AIC = NA, BIC = NA, Species = NA)

predictions_centromereModels$pred_centromere_arm[which(predictions_centromereModels$pred_centromere_arm == 
    Inf)] = NA

mod = lm(pred_telomere ~ centromeric_index, data = predictions_centromereModels)
summary(mod)$r.squared
```

```
## [1] 0.216888
```

```r
centromere_modselection$AdjustedRsquared[1] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[1] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[1] = AIC(mod)
centromere_modselection$BIC[1] = BIC(mod)

mod = lm(pred_centromere_arm ~ centromeric_index, data = predictions_centromereModels)
summary(mod)
```

```
## 
## Call:
## lm(formula = pred_centromere_arm ~ centromeric_index, data = predictions_centromereModels)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -127.606   -0.223    0.208    0.722   53.840 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)
## (Intercept)         0.8596     1.5918   0.540    0.589
## centromeric_index  -1.1815     3.2745  -0.361    0.718
## 
## Residual standard error: 9.299 on 422 degrees of freedom
##   (1 observation deleted due to missingness)
## Multiple R-squared:  0.0003084,	Adjusted R-squared:  -0.002061 
## F-statistic: 0.1302 on 1 and 422 DF,  p-value: 0.7184
```

```r
centromere_modselection$AdjustedRsquared[2] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[2] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[2] = AIC(mod)
centromere_modselection$BIC[2] = BIC(mod)

mod = lm(pred_centromere_chromosome ~ centromeric_index, data = predictions_centromereModels)
summary(mod)
```

```
## 
## Call:
## lm(formula = pred_centromere_chromosome ~ centromeric_index, 
##     data = predictions_centromereModels)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.44215 -0.06600  0.00659  0.08467  0.42330 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -0.01549    0.02349   -0.66     0.51    
## centromeric_index  1.00752    0.04826   20.88   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1375 on 423 degrees of freedom
## Multiple R-squared:  0.5075,	Adjusted R-squared:  0.5064 
## F-statistic: 435.9 on 1 and 423 DF,  p-value: < 2.2e-16
```

```r
centromere_modselection$AdjustedRsquared[3] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[3] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[3] = AIC(mod)
centromere_modselection$BIC[3] = BIC(mod)
```


Linear regression and model selection based on R2 criterion within species with at least 5 chromosomes.




```r
table(species$bestmodel)
```

```
## 
##  2  3 
##  7 30
```

```r
centromere_modselection$Species[1] = sum(species$bestmodel == 1)
centromere_modselection$Species[2] = sum(species$bestmodel == 2)
centromere_modselection$Species[3] = sum(species$bestmodel == 3)
```


Redo model selection but on a subset of chromosomes with no chromosome arm below 50cM.


```r
# Sample size: number of chromosomes:
nrow(predictions_centromereModels[which(predictions_centromereModels$pred_centromere_arm > 
    0), ])
```

```
## [1] 330
```



```r
predictions_centromereModels_subset = predictions_centromereModels[which(predictions_centromereModels$pred_centromere_arm > 
    0), ]

mod = lm(pred_telomere ~ centromeric_index, data = predictions_centromereModels_subset)
summary(mod)$r.squared
```

```
## [1] 0.1871755
```

```r
centromere_modselection$AdjustedRsquared[4] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[4] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[4] = AIC(mod)
centromere_modselection$BIC[4] = BIC(mod)

mod = lm(pred_centromere_arm ~ centromeric_index, data = predictions_centromereModels_subset)
summary(mod)
```

```
## 
## Call:
## lm(formula = pred_centromere_arm ~ centromeric_index, data = predictions_centromereModels_subset)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -1.762 -1.175 -0.848 -0.369 52.357 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)
## (Intercept)         0.8381     0.9167   0.914    0.361
## centromeric_index   1.5056     1.8734   0.804    0.422
## 
## Residual standard error: 4.541 on 328 degrees of freedom
## Multiple R-squared:  0.001965,	Adjusted R-squared:  -0.001077 
## F-statistic: 0.6459 on 1 and 328 DF,  p-value: 0.4222
```

```r
centromere_modselection$AdjustedRsquared[5] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[5] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[5] = AIC(mod)
centromere_modselection$BIC[5] = BIC(mod)

mod = lm(pred_centromere_chromosome ~ centromeric_index, data = predictions_centromereModels_subset)
summary(mod)
```

```
## 
## Call:
## lm(formula = pred_centromere_chromosome ~ centromeric_index, 
##     data = predictions_centromereModels_subset)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.43678 -0.06441  0.01070  0.07741  0.42520 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -0.009085   0.026649  -0.341    0.733    
## centromeric_index  0.988207   0.054461  18.145   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.132 on 328 degrees of freedom
## Multiple R-squared:  0.5009,	Adjusted R-squared:  0.4994 
## F-statistic: 329.2 on 1 and 328 DF,  p-value: < 2.2e-16
```

```r
centromere_modselection$AdjustedRsquared[6] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[6] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[6] = AIC(mod)
centromere_modselection$BIC[6] = BIC(mod)
```


Linear regression and model selection based on R2 criterion within species with at least 5 chromosomes.




```r
table(species_subset$bestmodel)
```

```
## 
##  2  3 
## 10 26
```

```r
centromere_modselection$Species[4] = sum(species_subset$bestmodel == 1)
centromere_modselection$Species[5] = sum(species_subset$bestmodel == 2)
centromere_modselection$Species[6] = sum(species_subset$bestmodel == 3)
```








Redo model selection but on a subset of chromosomes with a correct CI orientation assessed (remove wrong orientation) 50cM.


```r
id = which(paste(predictions_centromereModels$set, predictions_centromereModels$chromosome) %in% 
    paste(chromosome.stats$set[idx_wrongCIorientation], chromosome.stats$chromosome[idx_wrongCIorientation]))

paste(predictions_centromereModels$set[id], predictions_centromereModels$chromosome[id])
```

```
##  [1] "Boechera_stricta_Lee2020 2"                
##  [2] "Boechera_stricta_Lee2020 5"                
##  [3] "Brassica_napus_Yang2017 A07"               
##  [4] "Brassica_rapa_CodyMarkelz2017 A01"         
##  [5] "Capsella_rubella_Slotte2013 1"             
##  [6] "Capsella_rubella_Slotte2013 4"             
##  [7] "Capsicum_annuum_Han2016 3"                 
##  [8] "Capsicum_annuum_Han2016 12"                
##  [9] "Cenchrus_americanus_Pucher2017 2"          
## [10] "Cenchrus_americanus_Pucher2017 6"          
## [11] "Citrullus_lanatus_Ren2015 3"               
## [12] "Citrullus_lanatus_Ren2015 7"               
## [13] "Citrullus_lanatus_Ren2015 10"              
## [14] "Citrus_sinensis_Huang2018 4"               
## [15] "Coffea_canephora_Crouzillat2020 A"         
## [16] "Coffea_canephora_Crouzillat2020 K"         
## [17] "Cucumis_melo_Pereira2018 4"                
## [18] "Cucumis_melo_Pereira2018 7"                
## [19] "Cucumis_melo_Pereira2018 9"                
## [20] "Cucumis_sativus_Zhu2016 6"                 
## [21] "Cucumis_sativus_Zhu2016 5"                 
## [22] "Elaeis_guineensis_Yaakub2020 5"            
## [23] "Elaeis_guineensis_Yaakub2020 1"            
## [24] "Elaeis_guineensis_Yaakub2020 12"           
## [25] "Eucalyptus_grandis_Bertholome2015 1"       
## [26] "Eucalyptus_grandis_Bertholome2015 6"       
## [27] "Eucalyptus_grandis_Bertholome2015 7"       
## [28] "Eucalyptus_grandis_Bertholome2015 8"       
## [29] "Gossypium_hirsutum_Zhang2019 A11"          
## [30] "Lupinus_angustifolius_Zhou2017 8"          
## [31] "Lupinus_angustifolius_Zhou2017 10"         
## [32] "Lupinus_angustifolius_Zhou2017 19"         
## [33] "Lupinus_angustifolius_Zhou2017 15"         
## [34] "Lupinus_angustifolius_Zhou2017 3"          
## [35] "Lupinus_angustifolius_Zhou2017 7"          
## [36] "Malus_domestica_DiPierro2016 3"            
## [37] "Manihot_esculenta_ICGMC2015 4"             
## [38] "Manihot_esculenta_ICGMC2015 6"             
## [39] "Manihot_esculenta_ICGMC2015 8"             
## [40] "Manihot_esculenta_ICGMC2015 12"            
## [41] "Momordica_charantia_Mastumura2020 2"       
## [42] "Nelumbo_nucifera_Gui2018 1"                
## [43] "Nelumbo_nucifera_Gui2018 5"                
## [44] "Oryza_nivara_Ma2016 9"                     
## [45] "Oryza_nivara_Ma2016 12"                    
## [46] "Oryza_sativa_DeLeon2016 3"                 
## [47] "Oryza_sativa_DeLeon2016 5"                 
## [48] "Oryza_sativa_DeLeon2016 6"                 
## [49] "Phaseolus_vulgaris_Song2015 4"             
## [50] "Phaseolus_vulgaris_Song2015 10"            
## [51] "Phaseolus_vulgaris_Song2015 11"            
## [52] "Prunus_persica_Verde2017 3"                
## [53] "Sesamum_indicum_Wang2016 1"                
## [54] "Sesamum_indicum_Wang2016 3"                
## [55] "Sesamum_indicum_Wang2016 4"                
## [56] "Solanum_lycopersicum_Gonda2018 12"         
## [57] "Solanum_tuberosum_Endelman2016 1"          
## [58] "Solanum_tuberosum_Endelman2016 4"          
## [59] "Sorghum_bicolor_Zou2012 6"                 
## [60] "Sorghum_bicolor_Zou2012 7"                 
## [61] "Sorghum_bicolor_Zou2012 10"                
## [62] "Theobroma_cacao_Royaert2016 2"             
## [63] "Theobroma_cacao_Royaert2016 6"             
## [64] "Theobroma_cacao_Royaert2016 7"             
## [65] "Triticum_aestivum_GutierrezGonzalez2019 4B"
## [66] "Triticum_aestivum_GutierrezGonzalez2019 4D"
## [67] "Vitis_vinifera_Brault2020 4"               
## [68] "Vitis_vinifera_Brault2020 6"               
## [69] "Vitis_vinifera_Brault2020 9"               
## [70] "Vitis_vinifera_Brault2020 12"
```

```r
# Sample size: number of chromosomes:
nrow(predictions_centromereModels[-id, ])
```

```
## [1] 355
```



```r
predictions_centromereModels_subset = predictions_centromereModels[-id, ]

mod = lm(pred_telomere ~ centromeric_index, data = predictions_centromereModels_subset)
summary(mod)$r.squared
```

```
## [1] 0.3203657
```

```r
centromere_modselection$AdjustedRsquared[7] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[7] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[7] = AIC(mod)
centromere_modselection$BIC[7] = BIC(mod)

mod = lm(pred_centromere_arm ~ centromeric_index, data = predictions_centromereModels_subset)
summary(mod)
```

```
## 
## Call:
## lm(formula = pred_centromere_arm ~ centromeric_index, data = predictions_centromereModels_subset)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -127.484   -0.178    0.256    0.749   53.928 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)
## (Intercept)          0.907      1.915   0.474    0.636
## centromeric_index   -1.424      3.993  -0.357    0.722
## 
## Residual standard error: 10.16 on 352 degrees of freedom
##   (1 observation deleted due to missingness)
## Multiple R-squared:  0.000361,	Adjusted R-squared:  -0.002479 
## F-statistic: 0.1271 on 1 and 352 DF,  p-value: 0.7217
```

```r
centromere_modselection$AdjustedRsquared[8] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[8] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[8] = AIC(mod)
centromere_modselection$BIC[8] = BIC(mod)

mod = lm(pred_centromere_chromosome ~ centromeric_index, data = predictions_centromereModels_subset)
summary(mod)
```

```
## 
## Call:
## lm(formula = pred_centromere_chromosome ~ centromeric_index, 
##     data = predictions_centromereModels_subset)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.41421 -0.06353  0.00445  0.08306  0.38753 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       -0.05477    0.02459  -2.227   0.0266 *  
## centromeric_index  1.10408    0.05120  21.563   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1308 on 353 degrees of freedom
## Multiple R-squared:  0.5684,	Adjusted R-squared:  0.5672 
## F-statistic:   465 on 1 and 353 DF,  p-value: < 2.2e-16
```

```r
centromere_modselection$AdjustedRsquared[9] = summary(mod)$adj.r.squared
centromere_modselection$Pvalue[9] = summary(mod)$coefficients[2, 4]
centromere_modselection$AIC[9] = AIC(mod)
centromere_modselection$BIC[9] = BIC(mod)
```


Linear regression and model selection based on R2 criterion within species with at least 5 chromosomes.




```r
table(species_subset$bestmodel)
```

```
## 
##  2  3 
##  7 30
```

```r
centromere_modselection$Species[7] = sum(species_subset$bestmodel == 1)
centromere_modselection$Species[8] = sum(species_subset$bestmodel == 2)
centromere_modselection$Species[9] = sum(species_subset$bestmodel == 3)
```









Save model selection...





# Genomic landscapes correlate to recombination landscapes









```r
# Sample size - Species
sum(annotated$annotation)
```

```
## [1] 41
```

```r
# Sample size - Chromosomes
length(unique(data$map))
```

```
## [1] 483
```


```r
cat("Species with annotations that were integrated in the dataset:\n")
```

```
## Species with annotations that were integrated in the dataset:
```

```r
paste(annotated$species[annotated$annotation], annotated$accession[annotated$annotation])
```

```
##  [1] "Arabidopsis_thaliana GCA_000001735.2"            
##  [2] "Arachis_duranensis GCF_000817695.2"              
##  [3] "Arachis_hypogaea GCA_003713155.1"                
##  [4] "Brachypodium_distachyon GCA_000005505.4"         
##  [5] "Brassica_napus GCF_000686985.2"                  
##  [6] "Brassica_rapa GCF_000309985.1"                   
##  [7] "Camelina_sativa GCF_000633955.1"                 
##  [8] "Camellia_sinensis GCA_013676235.1"               
##  [9] "Capsella_rubella GCA_000375325.1"                
## [10] "Citrus_sinensis GCA_000317415.1"                 
## [11] "Cucumis_melo v4.0"                               
## [12] "Cucumis_sativus GCA_000004075.2"                 
## [13] "Cucurbita_maxima Cmaxima_v1.1"                   
## [14] "Cucurbita_pepo GCF_002806865.2"                  
## [15] "Dioscorea_alata GCA_002240015.2"                 
## [16] "Elaeis_guineensis GCF_000442705.1"               
## [17] "Eucalyptus_grandis Egrandis_297_v2.0"            
## [18] "Glycine_max GCF_000004515.4"                     
## [19] "Gossypium_hirsutum HAU_G.hirsutum_AD1genome_v1.1"
## [20] "Gossypium_raimondii GCA_000327365.1"             
## [21] "Helianthus_annuus GCA_002127325.1"               
## [22] "Hordeum_vulgare GCA_900075435.2"                 
## [23] "Lupinus_albus GCA_009771035.1"                   
## [24] "Malus_domestica GCA_004115385.1"                 
## [25] "Manihot_esculenta GCA_001659605.1"               
## [26] "Oryza_nivara GCA_000576065.1"                    
## [27] "Oryza_sativa GCA_001433935.1"                    
## [28] "Panicum_hallii GCA_002211085.2"                  
## [29] "Phaseolus_vulgaris GCA_001517995.1"              
## [30] "Prunus_mume GCF_000346735.1"                     
## [31] "Prunus_persica GCA_000346465.2"                  
## [32] "Sesamum_indicum GCF_000512975.1"                 
## [33] "Setaria_italica GCA_000263155.2"                 
## [34] "Solanum_lycopersicum GCA_000188115.3"            
## [35] "Solanum_tuberosum PGSC_DM_v4.03"                 
## [36] "Sorghum_bicolor GCA_000003195.3"                 
## [37] "Theobroma_cacao GCA_000403535.1"                 
## [38] "Triticum_aestivum GCA_900519105.1"               
## [39] "Vigna_unguiculata GCF_004118075.1"               
## [40] "Vitis_vinifera GCA_000003745.2"                  
## [41] "Zea_mays GCA_000005005.6"
```


```r
cat("Species without annotations:\n")
```

```
## Species without annotations:
```

```r
paste(annotated$species[!(annotated$annotation)], annotated$accession[!(annotated$annotation)])
```

```
##  [1] "Aegilops_speltoides GCA_000347335.2"  
##  [2] "Aegilops_speltoides GCA_000347335.2"  
##  [3] "Boechera_stricta NA"                  
##  [4] "Capsicum_annuum GCA_002878395.2"      
##  [5] "Cenchrus_americanus GCA_002174835.2"  
##  [6] "Citrullus_lanatus GCA_000238415.2"    
##  [7] "Coffea_canephora NA"                  
##  [8] "Draba_nivalis PRJNA657155"            
##  [9] "Juglans_regia NA"                     
## [10] "Lupinus_angustifolius GCA_002285895.2"
## [11] "Mangifera_indica GCA_011075055.1"     
## [12] "Momordica_charantia GCA_013281855.1"  
## [13] "Nelumbo_nucifera GCA_003033695.1"     
## [14] "Quercus_sp GCA_900291515.1"           
## [15] "Raphanus_sativus GCA_002197605.1"     
## [16] "Triticum_dicoccoides GCA_002575655.1" 
## [17] "Triticum_urartu GCA_003073215.1"
```

Most species had rarely more than 20 genes in a single window and variance of the averaged recombination rate  increased in the upper range of gene counts (we also observed an inflection of the variance around 20 genes in raw data), hence for gene count the quadratic regression was fitted to a subset of windows with at most 20 genes. We considered that windows with more genes were less reliable than the rest of the genome (e.g. too much annotations, overlapping genes).



```r
hist(data$rec.rate, breaks = 100, main = "", xlab = "Recombination rate", xlim = c(0, 
    20))
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-147-1.png)<!-- -->

```r
hist(data$gene_count, breaks = 100, main = "", xlab = "Gene count", xlim = c(0, 50))
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-147-2.png)<!-- -->

```r
cat("Proportion of windows with a null recombination rate (i.e. < 0.001):\n")
```

```
## Proportion of windows with a null recombination rate (i.e. < 0.001):
```

```r
sum(data$rec.rate < 0.001, na.rm = TRUE)/sum(!is.na(data$rec.rate), na.rm = TRUE)
```

```
## [1] 0.1259402
```




```r
ggplot(data = meanrecrate_genecount, aes(x = gene_count, group = species)) + geom_density(alpha = 0.4) + 
    labs(x = "Gene count") + geom_vline(aes(xintercept = 30)) + theme(axis.line = element_line(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
    panel.background = element_blank(), plot.title = element_text(color = "black", 
        size = 16, face = "bold.italic", hjust = 0.5), plot.subtitle = element_text(color = "black", 
        size = 16, hjust = 0.5), axis.title.x = element_text(color = "black", size = 16), 
    axis.title.y = element_text(color = "black", size = 16), axis.text = element_text(size = 16, 
        colour = "black"), strip.text = element_text(size = 12, colour = "black", 
        angle = 90), legend.position = "none", legend.key = element_rect(fill = "white", 
        size = 1), legend.key.width = unit(2, "cm"), legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16))
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-148-1.png)<!-- -->

## Within species, chromosome-level

At a chromosome level, we expect that windows with more genes will exhibit higher recombination rates. We tested this prediction with a Spearman correlation on 100kb windows for each chromosome. Furthermore, I tested the significance of the correlation for each chromosome and the significance of the relationship among species with a Linear Mixed Model (not shown here) taking into account repeated measurement within species and the structure of the dataset (species and chromosomes were random effects).


```
## [1] 0.4619748
```

```
##      2.5% 
## 0.4364954
```

```
##     97.5% 
## 0.4862635
```

<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-149-1.png" style="display: block; margin: auto;" />


```r
cat("Number of chromosomes analyzed:\n")
```

```
## Number of chromosomes analyzed:
```

```r
sum(!is.na(cor_genecount$correlation))
```

```
## [1] 447
```

```r
cat("Proportion of chromosomes with a significant Spearman correlation:\n")
```

```
## Proportion of chromosomes with a significant Spearman correlation:
```

```r
sum(cor_genecount$signif, na.rm = TRUE)/sum(!is.na(cor_genecount$signif))
```

```
## [1] 0.9082774
```




## Across species: Chromosomes pooled per species.

Is it a common pattern among species?


In order to assess the interspecific differences in the relationship between recombination rates and gene count, I pooled chromosomes per species.




```r
colnames(spcorr)[4] = "n"
kable(spcorr, digits = 2, align = "c", caption = "Species Spearman's rank correlation between recombination rates and gene number.")
```

<table>
<caption>(\#tab:unnamed-chunk-152)Species Spearman's rank correlation between recombination rates and gene number.</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> species </th>
   <th style="text-align:center;"> correlation </th>
   <th style="text-align:center;"> signif </th>
   <th style="text-align:center;"> n </th>
   <th style="text-align:center;"> phys.map.length.y </th>
   <th style="text-align:center;"> mean.recrate </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> 0.28 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 23.83 </td>
   <td style="text-align:center;"> 4.08 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> 0.53 </td>
   <td style="text-align:center;"> 0.88 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 100.95 </td>
   <td style="text-align:center;"> 1.02 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> 0.63 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 124.68 </td>
   <td style="text-align:center;"> 1.53 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> 0.62 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 45.62 </td>
   <td style="text-align:center;"> 6.05 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> 0.21 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 44.72 </td>
   <td style="text-align:center;"> 3.81 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> 0.34 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 26.14 </td>
   <td style="text-align:center;"> 4.59 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> 0.54 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 30.63 </td>
   <td style="text-align:center;"> 3.64 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> 0.89 </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 222.11 </td>
   <td style="text-align:center;"> 0.61 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> -0.07 </td>
   <td style="text-align:center;"> 0.57 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 15.84 </td>
   <td style="text-align:center;"> 3.15 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> 0.19 </td>
   <td style="text-align:center;"> 0.88 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 26.02 </td>
   <td style="text-align:center;"> 2.83 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> 0.58 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 30.99 </td>
   <td style="text-align:center;"> 4.27 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> 0.19 </td>
   <td style="text-align:center;"> 0.71 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 27.41 </td>
   <td style="text-align:center;"> 5.69 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> 0.27 </td>
   <td style="text-align:center;"> 0.58 </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 10.35 </td>
   <td style="text-align:center;"> 17.57 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> 0.27 </td>
   <td style="text-align:center;"> 0.50 </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 10.79 </td>
   <td style="text-align:center;"> 13.53 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> 0.11 </td>
   <td style="text-align:center;"> 0.67 </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 23.85 </td>
   <td style="text-align:center;"> 5.28 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> 0.83 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 46.94 </td>
   <td style="text-align:center;"> 2.65 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> 0.41 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 55.69 </td>
   <td style="text-align:center;"> 1.40 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> 0.74 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 48.13 </td>
   <td style="text-align:center;"> 3.11 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> 0.54 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 26 </td>
   <td style="text-align:center;"> 85.92 </td>
   <td style="text-align:center;"> 2.45 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> 0.62 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 57.63 </td>
   <td style="text-align:center;"> 3.00 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> 0.13 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 169.70 </td>
   <td style="text-align:center;"> 0.47 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> 0.40 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 654.86 </td>
   <td style="text-align:center;"> 0.24 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> 0.63 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 25 </td>
   <td style="text-align:center;"> 17.36 </td>
   <td style="text-align:center;"> 7.00 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> 0.64 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 38.38 </td>
   <td style="text-align:center;"> 2.00 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> 0.11 </td>
   <td style="text-align:center;"> 0.61 </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 28.80 </td>
   <td style="text-align:center;"> 6.12 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> 0.24 </td>
   <td style="text-align:center;"> 0.86 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 28.16 </td>
   <td style="text-align:center;"> 2.81 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> 0.51 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 31.10 </td>
   <td style="text-align:center;"> 4.42 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> 0.70 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 56.38 </td>
   <td style="text-align:center;"> 2.05 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> 0.76 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 46.80 </td>
   <td style="text-align:center;"> 2.00 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> 0.36 </td>
   <td style="text-align:center;"> 0.86 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 24.99 </td>
   <td style="text-align:center;"> 7.26 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> 0.71 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 28.21 </td>
   <td style="text-align:center;"> 2.85 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> 0.10 </td>
   <td style="text-align:center;"> 0.60 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 16.41 </td>
   <td style="text-align:center;"> 4.39 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> 0.67 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 44.59 </td>
   <td style="text-align:center;"> 3.61 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> 0.68 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 67.27 </td>
   <td style="text-align:center;"> 1.94 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> 0.52 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 60.42 </td>
   <td style="text-align:center;"> 1.27 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> 0.71 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 68.36 </td>
   <td style="text-align:center;"> 2.44 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> 0.83 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 33.05 </td>
   <td style="text-align:center;"> 2.60 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> 0.32 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 20 </td>
   <td style="text-align:center;"> 667.83 </td>
   <td style="text-align:center;"> 0.21 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> 0.72 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 43.04 </td>
   <td style="text-align:center;"> 1.78 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> 0.57 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 19 </td>
   <td style="text-align:center;"> 22.43 </td>
   <td style="text-align:center;"> 2.79 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> 0.47 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 210.63 </td>
   <td style="text-align:center;"> 0.87 </td>
  </tr>
</tbody>
</table>

```r
hist(spcorr$correlation, breaks = 30, xlab = "Species Spearman correlation", main = "", 
    xlim = c(-1, 1))
abline(v = 0, col = "Red", lty = 2)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-152-1.png)<!-- -->




```r
# df = spcorr[,1:4] colnames(df) = c('Species', 'Mean correlation', 'Proportion
# of significant correlations', 'Number of chromosomes') write.xlsx(x = df, file
# = paste(wd, 'tables/article_one_supp/tables_supplementary.xls', sep = ''),
# sheetName = 'S10_correlationRecombinationGenes', row.names = FALSE, append =
# TRUE)
```





Correlations seem uniformly distributed among species from zero to one. There is only few cases of negative correlation. The correlation previously inferred seems robust, yet the strength of the relationship vary greatly among species.
 
All chromosomes pooled, the mean recombination rate per gene count and confidence intervals at 95% were estimated for each species (1,000 bootstraps).



<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-154-1.png" alt="(ref:genecount-allspecies)"  />
<p class="caption">(\#fig:unnamed-chunk-154)(ref:genecount-allspecies)</p>
</div>

(ref:genecount-allspecies) (ref:genecount-allspecies) Recombination rate (cM/Mb) depends on the gene count for 44 species (chromosomes pooled per species). Mean recombination rates and confidence intervals (95%) estimated by 1,000 bootstraps. Gene count was made by counting the number of gene starting positions within the window.


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-155-1.png" alt="(ref:genecount-allspecies)"  />
<p class="caption">(\#fig:unnamed-chunk-155)(ref:genecount-allspecies)</p>
</div>

We see that the strength of the pattern varies among species. Intuitively, we can think that the slope of the relationship may depend on chromosome size, as many other patterns already observed.


```r
# cor_genecount = merge(cor_genecount, chromosome.stats[,c(1, 3, 8)]) plot(x =
# cor_genecount$phys.map.length, y = cor_genecount$correlation, log = 'x')

# Estimate the slope of the relationship for each species Data = 100kb windows

# Estimate the slope of the relationship for each species Data = pooled
# chromosomes, recombination rate per gene count
spcorr$intercept = NA
spcorr$slope = NA
for (i in 1:nrow(spcorr)) {
    subset = meanrecrate_genecount[which(gsub("_", " ", meanrecrate_genecount$species) == 
        spcorr$species[i]), ]
    lm_model = lm(mean_rec ~ gene_count, data = subset)
    spcorr$intercept[i] = lm_model$coefficients[1]
    spcorr$slope[i] = lm_model$coefficients[2]
}

# predictions not verified None effect of chromosome length
plot(x = spcorr$phys.map.length, y = spcorr$slope, log = "x")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-156-1.png)<!-- -->

```r
cor.test(x = spcorr$phys.map.length, y = spcorr$slope, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  spcorr$phys.map.length and spcorr$slope
## S = 11346, p-value = 0.9424
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## 0.01167247
```

```r
# And mean recombination rate?
plot(x = spcorr$mean.recrate, y = spcorr$slope, log = "x", xlab = "Mean recombination rate", 
    ylab = "Slope")
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-156-2.png)<!-- -->

```r
# Significance of teh relatiionship
cor.test(x = spcorr$mean.recrate, y = spcorr$slope, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  spcorr$mean.recrate and spcorr$slope
## S = 10572, p-value = 0.6219
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## 0.07909408
```

```r
summary(lm(slope ~ mean.recrate, data = spcorr))
```

```
## 
## Call:
## lm(formula = slope ~ mean.recrate, data = spcorr)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.25164 -0.08751 -0.03233  0.05596  0.41274 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  0.099649   0.032511   3.065  0.00394 **
## mean.recrate 0.001453   0.006660   0.218  0.82839   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1385 on 39 degrees of freedom
## Multiple R-squared:  0.00122,	Adjusted R-squared:  -0.02439 
## F-statistic: 0.04762 on 1 and 39 DF,  p-value: 0.8284
```


## Is it a common pattern among species?






```r
summary(lm(mean_rec ~ gene_count, data = meanrecrate_genecount_subset))
```

```
## 
## Call:
## lm(formula = mean_rec ~ gene_count, data = meanrecrate_genecount_subset)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -3.4534 -0.3457  0.0101  0.3602  3.1616 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -1.219824   0.042469  -28.72   <2e-16 ***
## gene_count   0.124219   0.003695   33.62   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.6328 on 820 degrees of freedom
##   (4 observations deleted due to missingness)
## Multiple R-squared:  0.5795,	Adjusted R-squared:  0.579 
## F-statistic:  1130 on 1 and 820 DF,  p-value: < 2.2e-16
```

Quadratic regression.


```r
summary(lm(mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount_subset))
```

```
## 
## Call:
## lm(formula = mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount_subset)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -3.1784 -0.3209 -0.0092  0.3239  3.5239 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     -1.5821115  0.0572722 -27.624   <2e-16 ***
## gene_count       0.2394674  0.0133397  17.952   <2e-16 ***
## I(gene_count^2) -0.0058240  0.0006501  -8.959   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.6043 on 819 degrees of freedom
##   (4 observations deleted due to missingness)
## Multiple R-squared:  0.6171,	Adjusted R-squared:  0.6161 
## F-statistic: 659.9 on 2 and 819 DF,  p-value: < 2.2e-16
```



<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-160-1.png" alt="(ref:genecount-quadraticregression)"  />
<p class="caption">(\#fig:unnamed-chunk-160)(ref:genecount-quadraticregression)</p>
</div>

(ref:genecount-quadraticregression) (ref:genecount-quadraticregression) Recombination rate depends on the gene count for 44 species both across and within species (A) Mean recombination rates (chromosomes pooled per species, one color per species) estimated by 1,000 bootstraps and standardized within species. Gene count was made by counting the number of gene starting positions within 100kb windows with a recombination rate. Black line with grey ribbon is the quadratic regression estimated by linear regression. (B) Distribution of chromosome Spearman rank correlations between the number of genes and the recombination rate in 100kb windows (n = 488 chromosomes). The black vertical lines are the mean correlation across chromosomes (solid) and 95% confidence interval (dashed) estimated by 1,000 bootstraps. (C) Slope of the species linear regression between gene count and recombination rates. Linear regression estimated on chromosome-level data with recombination rates estimated in quantiles of gene count. The linear regression is not significant (p = 0.72).


#### Figure 8







Though the relationship seems a common feature among species, some of them are outliers and show a stronger relationship than the rest of species. It would be interesting to look at specific relationships to identify outlying patterns.



## Gene landscapes correlate to recombination landscape

Plot gene landscapes in the same manner as recombination landscapes in (ref:DistancesRelative-pooledChromosomes-patterns).

Number of genes ~ distance to the telomere.







<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-163-1.png" alt="(ref:DistancesRelative-pooledChromosomes-GeneCount)"  />
<p class="caption">(\#fig:unnamed-chunk-163)(ref:DistancesRelative-pooledChromosomes-GeneCount)</p>
</div>

(ref:DistancesRelative-pooledChromosomes-GeneCount) Standardized gene count as a function of the relative distance from the telomere along the chromosome (physical distances expressed in 20 bins). Chromosomes were split in halves, distance 0.5 being the center of the chromosome. Then, chromosomes were pooled per species. Each color is a species. Linear quadratic regression estimated is the black line (95% confidence interval in dark grey).



<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-164-1.png" alt="(ref:DistancesRelative-pooledChromosomes-LOESS-GeneCount)"  />
<p class="caption">(\#fig:unnamed-chunk-164)(ref:DistancesRelative-pooledChromosomes-LOESS-GeneCount)</p>
</div>

(ref:DistancesRelative-pooledChromosomes-LOESS-GeneCount) Standardized gene count as a function of the relative distance from the telomere along the chromosome (physical distances expressed in 20 bins). Chromosomes were split in halves, distance 0.5 being the center of the chromosome. Then, chromosomes were pooled per species. Each color is a species. A loess regression was estimated for each species.


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-165-1.png" alt="(ref:DistancesRelative-pooledChromosomes-genecount)"  />
<p class="caption">(\#fig:unnamed-chunk-165)(ref:DistancesRelative-pooledChromosomes-genecount)</p>
</div>

(ref:DistancesRelative-pooledChromosomes-genecount) (ref:DistancesRelative-pooledChromosomes-genecount)


#### Figure S7





Classify gene patterns as we did for crossover patterns.







```r
table(telomere_pattern$gene_pattern)
```

```
## 
## exception      peak  telomere 
##         5         7        29
```

```r
table(telomere_pattern$pattern)
```

```
## 
## exception      peak  telomere 
##         7        16        34
```

```r
# Species that changed of pattern among recombination/genes
telomere_pattern$species[which(telomere_pattern$pattern != telomere_pattern$gene_pattern)]
```

```
##  [1] "Brassica napus"     "Dioscorea alata"    "Eucalyptus grandis"
##  [4] "Malus domestica"    "Mangifera indica"   "Manihot esculenta" 
##  [7] "Panicum hallii"     "Prunus persica"     "Sesamum indicum"   
## [10] "Sorghum bicolor"    "Vitis vinifera"
```



```r
knitr::kable(telomere_pattern, label = c("Species", "Recombination pattern", "Gene density pattern"), 
    caption = "genelandscapes-patterns", align = "c")
```

<table>
<caption>(\#tab:SpeciesRecombination patternGene density pattern)genelandscapes-patterns</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> species </th>
   <th style="text-align:center;"> pattern </th>
   <th style="text-align:center;"> gene_pattern </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Boechera stricta </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsicum annuum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cenchrus americanus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrullus lanatus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Coffea canephora </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Draba nivalis </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Juglans regia </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Mangifera indica </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Momordica charantia </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Nelumbo nucifera </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Quercus sp </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Raphanus sativus </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> exception </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum dicoccoides </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum urartu </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> peak </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
  </tr>
</tbody>
</table>




<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-171-1.png" alt="(ref:DistancesRelative-pooledChromosomes-patterns-genecount)"  />
<p class="caption">(\#fig:unnamed-chunk-171)(ref:DistancesRelative-pooledChromosomes-patterns-genecount)</p>
</div>

(ref:DistancesRelative-pooledChromosomes-patterns-genecount) (ref:DistancesRelative-pooledChromosomes-patterns-genecount) Gene counts are higher in distal regions and lower near the center of the chromosome. Standardized gene count as a function of the relative distance from the telomere along the chromosome (physical distances expressed in 20 bins). Two patterns were identified based on the distribution of recombination rates along the telomere-centromere axis and species were pooled accordingly. Recombination rates decrease immediately from the tip of the chromosome for 35 species (dark grey line and ribbon) or recombination rate is reduced in telomeric regions and the peak of recombination is in a sub-telomeric region for 16 species (light grey line and ribbon). The solid line represents the mean gene count estimated in a bin and upper and lower boundaries of the ribbon represent the maximum and minimum values for a particular pattern. Besides, patterns that were not classified (6 species) were represented by loess regression in grey dashed lines. For estimating gene counts in bins of relative distance, chromosomes were split in halves, a distance of 0.5 being the center of the chromosome. Then, chromosomes were pooled per species.

#### Figure 9





# Genetic shuffling



## Correlations with genetic shuffling rates


```echo
# Sample size
length(!is.na(chromosome.stats.shuffling$geneticshuffling_genedist))
```


```r
hist(sqrt(chromosome.stats.shuffling$peripherybias_ratio), main = "", xlab = "Periphery-bias ratio (square root)", 
    breaks = 40)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-174-1.png)<!-- -->

```r
hist(chromosome.stats.shuffling$geneticshuffling, main = "", xlab = "Genetic shuffling", 
    breaks = 40)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-174-2.png)<!-- -->

```r
hist(chromosome.stats.shuffling$geneticshuffling_genedist, main = "", xlab = "Genetic shuffling", 
    breaks = 40)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-174-3.png)<!-- -->


### Chromosomal level

Hypothesis: Crossovers clustered in distal regions are less efficient than crossovers evenly distributed in the chromosome.


```r
lmer.model = lmer(geneticshuffling ~ (peripherybias_ratio) + (1 | species), data = chromosome.stats.shuffling)
summary(lmer.model)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: geneticshuffling ~ (peripherybias_ratio) + (1 | species)
##    Data: chromosome.stats.shuffling
## 
## REML criterion at convergence: -1167.8
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.1204 -0.5178 -0.0442  0.4157  6.8858 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  species  (Intercept) 0.005632 0.07505 
##  Residual             0.002656 0.05153 
## Number of obs: 421, groups:  species, 40
## 
## Fixed effects:
##                       Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)           0.255838   0.013177  49.996897  19.416  < 2e-16 ***
## peripherybias_ratio  -0.013249   0.002194 408.850703  -6.038  3.5e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## prphrybs_rt -0.362
```

```r
r.squaredGLMM(lmer.model)
```

```
##             R2m       R2c
## [1,] 0.05600436 0.6975317
```

```r
# Confidence interval of the coefficients was computed by bootstrap Parametric
# bootstrap
lmer.confint = confint(lmer.model, level = 0.95, method = "boot", nsim = 1000, boot.type = "norm")
```

```
## Computing bootstrap confidence intervals ...
```

```r
# Quadratic function and C.I.
lmerfun = function(x, c) {
    fixef(lmer.model)[[1]] + fixef(lmer.model)[[2]] * x
}
lowerlmerfun = function(x, c) {
    lmer.confint[[3]] + lmer.confint[[4]] * x
}
upperlmerfun = function(x, c) {
    lmer.confint[[7]] + lmer.confint[[8]] * x
}
# The C.I. values
df.new = data.frame(peripherybias_ratio = (chromosome.stats.shuffling$peripherybias_ratio))
df.new$lwr.pred = lowerlmerfun((chromosome.stats.shuffling$peripherybias_ratio))
df.new$upr.pred = upperlmerfun((chromosome.stats.shuffling$peripherybias_ratio))
```
<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-176-1.png" alt="(ref:shuffling-peripherybias)"  />
<p class="caption">(\#fig:unnamed-chunk-176)(ref:shuffling-peripherybias)</p>
</div>


(ref:shuffling-peripherybias) (ref:shuffling-peripherybias)

#### Figure S9




As expected intuitively, larger genetic maps have a higher Genetic shuffling.



```r
lmer.model = lmer(geneticshuffling ~ linkage.map.length.correctedHW + (1 | species), 
    data = chromosome.stats.shuffling)
summary(lmer.model)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: geneticshuffling ~ linkage.map.length.correctedHW + (1 | species)
##    Data: chromosome.stats.shuffling
## 
## REML criterion at convergence: -1488.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.0971 -0.4578 -0.0038  0.3746  7.0112 
## 
## Random effects:
##  Groups   Name        Variance Std.Dev.
##  species  (Intercept) 0.004415 0.06644 
##  Residual             0.001260 0.03550 
## Number of obs: 430, groups:  species, 40
## 
## Fixed effects:
##                                 Estimate Std. Error        df t value Pr(>|t|)
## (Intercept)                    5.443e-02  1.319e-02 8.042e+01   4.127 8.93e-05
## linkage.map.length.correctedHW 1.393e-03  6.192e-05 4.238e+02  22.489  < 2e-16
##                                   
## (Intercept)                    ***
## linkage.map.length.correctedHW ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## lnkg.mp..HW -0.582
```

```r
r.squaredGLMM(lmer.model)
```

```
##            R2m       R2c
## [1,] 0.4294228 0.8732853
```


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-179-1.png" alt="(ref:shuffling-linkagemaplength)"  />
<p class="caption">(\#fig:unnamed-chunk-179)(ref:shuffling-linkagemaplength)</p>
</div>
(ref:shuffling-linkagemaplength) (ref:shuffling-linkagemaplength)

#### Figure S8






### Species level

We need information about chromosome relative sizes and all chromosomes to compute the chromosomal component Genetic shuffling for a species.



## Gene distances instead of physical distances in Mb

#### Figure S10

Save gene densities for information in a Supplementary pdf...








#### Figure S11

Save Marey maps with gene distances in a Supplementary pdf...








Some examples of well chosen recombination landscapes in physical/gene distance.


<!-- ```{r echo = FALSE} -->
<!-- # Produce one figure per Marey map with gene distances instead of physical distances in Mb. -->
<!-- # List of all maps (all maps with a gc_genes file giving the list of genes and their position in a reference genome) -->
<!-- list = system(paste("ls ", wd, "/data-cleaned/genome/gene_distance/", sep = ""), intern = TRUE) -->
<!-- list = list[grep("*.txt.gz", list)] -->
<!-- list = gsub("_gene_positions.txt.gz", "", list) -->
<!-- # list -->
<!-- for (sp in list) { -->
<!--   # cat("================================\n") -->
<!--   # cat("Processing", sp, "...\n") -->
<!--   # Load the map with gene distances and genetic distances (cM) -->
<!--   gene_distance_data = read.table(file = gzfile(paste(wd, "data-cleaned/genome/gene_distance/", sp, sep = "")), header = TRUE, sep = "\t", stringsAsFactors = FALSE) -->
<!--   # Run for each chromosome -->
<!--   list_chr = unique(gene_distance_data$map, na.rm = TRUE) -->
<!--   list_chr = list_chr[which(!is.na(list_chr))] -->
<!--   for (chr in list_chr) { -->
<!--     gene_distance_map = gene_distance_data[gene_distance_data$map == chr,] -->
<!--     # Some maps have gene distances == 0 -->
<!--     if (max(gene_distance_map$gene_distance, na.rm = TRUE) > 0) { -->
<!--       # Save the figure -->
<!--       set = unique(gene_distance_map$set) -->
<!--       set = set[which(!is.na(set))] -->
<!--       # cat("Saving", set, "chromosome", chr, "\n") -->
<!--     png(filename = paste(wd, "output/gene_distances/marey_maps/", set, "_chr", chr, ".png", sep = ""), width = 1600, height = 860) -->
<!--     par(mfrow = c(1,2)) -->
<!--     plot(gene_distance_map$phys/1000000, gene_distance_map$gen, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)", main = paste(set, "_chr", chr, sep = "")) -->
<!--     lines(x = c(min(gene_distance_map$phys/1000000, na.rm = TRUE), max(gene_distance_map$phys/1000000, na.rm = TRUE)), y = c(min(gene_distance_map$gen, na.rm = TRUE), max(gene_distance_map$gen, na.rm = TRUE))) -->
<!--     plot(gene_distance_map$gene_distance, gene_distance_map$gen, xlab = "Gene distance\n(cumulative number of genes)", ylab = "Genetic distance (cM)", main = paste(set, "_chr", chr, sep = "")) -->
<!--     lines(x = c(min(gene_distance_map$gene_distance, na.rm = TRUE), max(gene_distance_map$gene_distance, na.rm = TRUE)), y = c(min(gene_distance_map$gen, na.rm = TRUE), max(gene_distance_map$gen, na.rm = TRUE))) -->
<!--     par(mfrow = c(1,1)) -->
<!--     dev.off() -->
<!--     # Interestingly, -->
<!--     png(filename = paste(wd, "output/gene_distances/map_correlations/", set, "_chr", chr, ".png", sep = ""), width = 860, height = 860) -->
<!--     plot(gene_distance_map$phys/1000000, gene_distance_map$gene_distance, xlab = "Physical distance (Mb)", ylab = "Gene distance", main = paste(set, "_chr", chr, sep = "")) -->
<!--     lines(x = c(0, max(gene_distance_map$phys/1000000, na.rm = TRUE)), y = c(0, max(gene_distance_map$gene_distance, na.rm = TRUE))) -->
<!--     dev.off() -->
<!--     } -->
<!--   } -->
<!-- } -->
<!-- ``` -->



<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-185-1.png" alt="(ref:marey-genedistances)"  />
<p class="caption">(\#fig:unnamed-chunk-185)(ref:marey-genedistances)</p>
</div>


(ref:marey-genedistances) (ref:marey-genedistances landscape_capsella, landscape_arabidopsis, landscape_malus, landscape_eucalyptus, landscape_oryza, landscape_zea)

#### Figure 10








Genetic shuffling should be (1) higher and (2) more homogeneous when estimated in gene distances rather than physical distances.


```r
# Sample size
nrow(chromosome.stats.shuffling)
```

```
## [1] 430
```

```r
length(unique(chromosome.stats.shuffling$species))
```

```
## [1] 40
```

```r
summary(chromosome.stats.shuffling$geneticshuffling)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0323  0.1602  0.2095  0.2200  0.2676  0.7157
```

```r
summary(chromosome.stats.shuffling$geneticshuffling_genedist)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1021  0.2055  0.2509  0.2655  0.2961  0.9186
```

```r
# The strength of the difference in genetic shuffling rate between bp and gene
# distances
mean(chromosome.stats.shuffling$geneticshuffling_genedist, na.rm = TRUE) - mean(chromosome.stats.shuffling$geneticshuffling, 
    na.rm = TRUE)
```

```
## [1] 0.04546639
```

```r
(mean(chromosome.stats.shuffling$geneticshuffling_genedist, na.rm = TRUE)/mean(chromosome.stats.shuffling$geneticshuffling, 
    na.rm = TRUE)) - 1
```

```
## [1] 0.2066655
```

```r
# ggplot(data = chromosome.stats.shuffling, aes(x = geneticshuffling)) +
# geom_histogram(color = 'Red', alpha = 0.5, bins = 50) + geom_histogram(aes(x =
# geneticshuffling_genedist), color = 'Blue', alpha = 0.5, bins = 50)

ggplot(data = chromosome.stats.shuffling, aes(x = geneticshuffling)) + geom_density(color = "Red", 
    alpha = 0.5) + geom_density(aes(x = geneticshuffling_genedist), color = "Blue", 
    alpha = 0.5)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-187-1.png)<!-- -->

```r
boxplot(chromosome.stats.shuffling$geneticshuffling, chromosome.stats.shuffling$geneticshuffling_genedist, 
    names = c("Physical distances", "Gene distances"))
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-187-2.png)<!-- -->

```r
plot(chromosome.stats.shuffling$geneticshuffling, chromosome.stats.shuffling$geneticshuffling_genedist, 
    xlab = "Genetic shuffling", ylab = "Genetic shuffling in gene distances")
abline(0, 1)
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-187-3.png)<!-- -->


<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-188-1.png" alt="(ref:dist-shufflingrates)"  />
<p class="caption">(\#fig:unnamed-chunk-188)(ref:dist-shufflingrates)</p>
</div>

I made the difference in genetic shuffling rate when changing genomic distances for gene distances (i.e. geneticshuffling_genedist - geneticshuffling). A positive value indicates that considering gene distances increases the genetic shuffling rate.



```r
# Permutation test of the significance of the difference in genetic shuffling
# rate
nboot = 1000
boot = numeric(nboot)
set.seed(42)
for (i in 1:nboot) {
    boot[i] = mean(sample((chromosome.stats.shuffling$geneticshuffling_genedist - 
        chromosome.stats.shuffling$geneticshuffling), replace = TRUE), na.rm = TRUE)
}
(mdiff = mean(boot))  # The mean difference in genetic shuffling rate (phys - gene distance)
```

```
## [1] 0.04540133
```

```r
# i.e. a value under 0 means a higher genetic shuffling rate when considering
# gene distances Is it significant at 0.05
pthreshold = 0.001
(mquant = quantile(boot, c(pthreshold/2, 1 - pthreshold/2)))
```

```
##      0.05%     99.95% 
## 0.03843896 0.05311997
```

```r
# The difference is small (though the original shuffling rate is not a large
# value) But the difference is highly significant
```




```
## 
## Call:
## lm(formula = diff ~ chrsize, data = df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.07294 -0.02894 -0.01314  0.01741  0.15455 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 3.274e-02  8.936e-03   3.663 0.000756 ***
## chrsize     1.792e-04  5.424e-05   3.303 0.002089 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.04847 on 38 degrees of freedom
## Multiple R-squared:  0.2231,	Adjusted R-squared:  0.2026 
## F-statistic: 10.91 on 1 and 38 DF,  p-value: 0.002089
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-189-1.png" alt="(ref:differences-shufflingrates)"  />
<p class="caption">(\#fig:unnamed-chunk-189)(ref:differences-shufflingrates)</p>
</div>
(ref:dist-shufflingrates) (ref:dist-shufflingrates)




```r
wilcox.test(chromosome.stats.shuffling$geneticshuffling, chromosome.stats.shuffling$geneticshuffling_genedist)
```

```
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  chromosome.stats.shuffling$geneticshuffling and chromosome.stats.shuffling$geneticshuffling_genedist
## W = 64418, p-value = 1.404e-14
## alternative hypothesis: true location shift is not equal to 0
```

Longer chromosomes have a higher difference between genetic shuffling rates in bp/gene distances.




```r
lmer.model = lmer((geneticshuffling_genedist - geneticshuffling) ~ log10(phys.map.length) + 
    (1 | species), data = chromosome.stats.shuffling)
summary(lmer.model)
```

```
## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
## lmerModLmerTest]
## Formula: 
## (geneticshuffling_genedist - geneticshuffling) ~ log10(phys.map.length) +  
##     (1 | species)
##    Data: chromosome.stats.shuffling
## 
## REML criterion at convergence: -2032.8
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.5352 -0.5316  0.0722  0.5623  3.6205 
## 
## Random effects:
##  Groups   Name        Variance  Std.Dev.
##  species  (Intercept) 0.0018671 0.04321 
##  Residual             0.0003497 0.01870 
## Number of obs: 430, groups:  species, 40
## 
## Fixed effects:
##                          Estimate Std. Error         df t value Pr(>|t|)    
## (Intercept)             -0.045069   0.016333 189.910554  -2.759  0.00636 ** 
## log10(phys.map.length)   0.055999   0.008915 298.406164   6.281 1.19e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr)
## lg10(phy..) -0.905
```

```r
r.squaredGLMM(lmer.model)
```

```
##           R2m       R2c
## [1,] 0.209811 0.8753577
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-192-1.png" alt="(ref:shuffling-differencesshufflingrates)"  />
<p class="caption">(\#fig:unnamed-chunk-192)(ref:shuffling-differencesshufflingrates)</p>
</div>


(ref:shuffling-differencesshufflingrates) (ref:shuffling-differencesshufflingrates)


## Convergence among patterns?



Is there convergence between patterns, i.e. between telomere patterns, gene patterns, gene distance maps and genetic shuffling?



```r
pattern_convergence = telomere_pattern
# For each species, are gene distance maps more/less/similarly homogeneous than
# genomic distances?  For each chromosome, compute the RMSE for gene distances
# and genomic distances RMSE = sqrt(sum((Xobs - Xpred)^2))
chromosome.stats.shuffling$RMSEgenomic = NA
chromosome.stats.shuffling$RMSEgene = NA
for (i in 1:nrow(chromosome.stats.shuffling)) {
    map = read.table(paste(wd, "data-cleaned/genome/gene_distance/", chromosome.stats.shuffling$set[i], 
        ".txt.gz", sep = ""), header = TRUE)
    map = map[which(map$map == chromosome.stats.shuffling$chromosome[i]), ]
    # Relative distances
    map$phys = map$phys/max(map$phys, na.rm = TRUE)
    map$gene_distance = map$gene_distance/max(map$gene_distance, na.rm = TRUE)
    # Estimate predicted (i.e. diagonal line) The formula for a strictly linear line
    # is y = 0 + (max(y)/1)*x Because x is a relative distance (O to 1) So slope
    # is...
    slope = max(map$gen, na.rm = TRUE)
    # And predicted...
    map$predictedGenomic = slope * map$phys
    map$predictedGene = slope * map$gene_distance
    # Compute RMSE
    chromosome.stats.shuffling$RMSEgenomic[i] = sqrt(sum((map$gen - map$predictedGenomic)^2))
    chromosome.stats.shuffling$RMSEgene[i] = sqrt(sum((map$gen - map$predictedGene)^2))
    rm(map)
}
chromosome.stats.shuffling$diffRMSE = chromosome.stats.shuffling$RMSEgene - chromosome.stats.shuffling$RMSEgenomic

# for each species, estimate the mean difference between RMSEs And classify as
# less/more homogeneous (>0/<0)
pattern_convergence$diffRMSE = NA
pattern_convergence$genedist_pattern = NA
for (i in 1:nrow(pattern_convergence)) {
    pattern_convergence$diffRMSE[i] = mean(chromosome.stats.shuffling$diffRMSE[which(chromosome.stats.shuffling$species == 
        pattern_convergence$species[i])], na.rm = TRUE)

    if (pattern_convergence$diffRMSE[i] == "NaN") {
        pattern_convergence$diffRMSE[i] = NA
    } else {
        if (pattern_convergence$diffRMSE[i] > 0) {
            pattern_convergence$genedist_pattern[i] = "less"
        } else {
            if (pattern_convergence$diffRMSE[i] < 0) {
                pattern_convergence$genedist_pattern[i] = "more"
            }
        }
    }

}


# For each species, are differences in genetic shuffling positive (increased) or
# negative (decreased)?  geneticshuffling_genedist - geneticshuffling
pattern_convergence$diffgeneticshuffl = NA
for (i in 1:nrow(pattern_convergence)) {
    pattern_convergence$diffgeneticshuffl[i] = mean(chromosome.stats.shuffling$geneticshuffling_genedist[which(chromosome.stats.shuffling$species == 
        pattern_convergence$species[i])] - chromosome.stats.shuffling$geneticshuffling[which(chromosome.stats.shuffling$species == 
        pattern_convergence$species[i])], na.rm = TRUE)
    if (pattern_convergence$diffgeneticshuffl[i] == "NaN") {
        pattern_convergence$diffgeneticshuffl[i] = NA
    }
}

# Mean chromosome size
pattern_convergence$meanchrsize = NA
for (i in 1:nrow(pattern_convergence)) {
    pattern_convergence$meanchrsize[i] = mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == 
        pattern_convergence$species[i])], na.rm = TRUE)
}
```



```r
# Number of maps classified as more homogeneous in gene distances:
sum(chromosome.stats.shuffling$diffRMSE < 0, na.rm = TRUE)
```

```
## [1] 303
```

```r
# Proportion:
sum(chromosome.stats.shuffling$diffRMSE < 0, na.rm = TRUE)/sum(!is.na(chromosome.stats.shuffling$diffRMSE), 
    na.rm = TRUE)
```

```
## [1] 0.7481481
```

```r
# Total:
sum(!is.na(chromosome.stats.shuffling$diffRMSE), na.rm = TRUE)
```

```
## [1] 405
```


```r
# Cross-tables for checking convergence
table(pattern_convergence$gene_pattern, pattern_convergence$pattern)
```

```
##            
##             exception peak telomere
##   exception         3    1        1
##   peak              0    5        2
##   telomere          1    6       22
```

```r
# How many species more and less?
table(pattern_convergence$genedist_pattern)
```

```
## 
## less more 
##    9   30
```

```r
# Interaction telomere pattern and gene distance
table(pattern_convergence$gene_pattern, pattern_convergence$genedist_pattern)
```

```
##            
##             less more
##   exception    3    2
##   peak         4    3
##   telomere     2   25
```

```r
# Species that changed of pattern among recombination/genes
pattern_convergence$species[which(pattern_convergence$pattern != pattern_convergence$gene_pattern)]
```

```
##  [1] "Brassica napus"     "Dioscorea alata"    "Eucalyptus grandis"
##  [4] "Malus domestica"    "Mangifera indica"   "Manihot esculenta" 
##  [7] "Panicum hallii"     "Prunus persica"     "Sesamum indicum"   
## [10] "Sorghum bicolor"    "Vitis vinifera"
```


```r
# Cross-tables for checking convergence
table(pattern_convergence$gene_pattern)
```

```
## 
## exception      peak  telomere 
##         5         7        29
```

```r
table(pattern_convergence$genedist_pattern)
```

```
## 
## less more 
##    9   30
```

```r
chisq.test(table(pattern_convergence$genedist_pattern, pattern_convergence$gene_pattern))
```

```
## Warning in chisq.test(table(pattern_convergence$genedist_pattern,
## pattern_convergence$gene_pattern)): Chi-squared approximation may be incorrect
```

```
## 
## 	Pearson's Chi-squared test
## 
## data:  table(pattern_convergence$genedist_pattern, pattern_convergence$gene_pattern)
## X-squared = 12.151, df = 2, p-value = 0.002299
```




```r
levels(pattern_convergence$pattern) = c("exception", "sub-telomere", "telomere")

knitr::kable(pattern_convergence, label = c("Species", "Recombination pattern", "Gene density pattern", 
    "RMSE difference", "Homogeneization effect of gene distance", "Genetic shuffling difference", 
    "Mean chromosome size"), caption = "genelandscapes-patterns", align = "c")
```

<table>
<caption>(\#tab:SpeciesRecombination patternGene density patternRMSE differenceHomogeneization effect of gene distanceGenetic shuffling differenceMean chromosome size)genelandscapes-patterns</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> species </th>
   <th style="text-align:center;"> pattern </th>
   <th style="text-align:center;"> gene_pattern </th>
   <th style="text-align:center;"> diffRMSE </th>
   <th style="text-align:center;"> genedist_pattern </th>
   <th style="text-align:center;"> diffgeneticshuffl </th>
   <th style="text-align:center;"> meanchrsize </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> Aegilops speltoides </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 517.44324 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arabidopsis thaliana </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> 30.5745792 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> 0.0064004 </td>
   <td style="text-align:center;"> 23.82927 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis duranensis </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 0.0579436 </td>
   <td style="text-align:center;"> 100.95097 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Arachis hypogaea </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -614.0625554 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.2072079 </td>
   <td style="text-align:center;"> 124.68497 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Boechera stricta </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 26.52555 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brachypodium distachyon </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -132.9239396 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0320681 </td>
   <td style="text-align:center;"> 45.62172 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica napus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> -24.6280031 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0224256 </td>
   <td style="text-align:center;"> 44.71564 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Brassica rapa </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> -64.6465361 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0371777 </td>
   <td style="text-align:center;"> 26.14066 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camelina sativa </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -21.9977716 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0098442 </td>
   <td style="text-align:center;"> 30.63062 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Camellia sinensis </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> -0.7685128 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0128694 </td>
   <td style="text-align:center;"> 222.10960 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsella rubella </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> 28.7411007 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> -0.0051357 </td>
   <td style="text-align:center;"> 15.83906 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Capsicum annuum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 227.12175 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cenchrus americanus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 223.50536 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrullus lanatus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 33.59728 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Citrus sinensis </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 26.02017 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Coffea canephora </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 53.00321 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis melo </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -76.2766136 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0580861 </td>
   <td style="text-align:center;"> 30.99121 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucumis sativus </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> 17.0914585 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> 0.0256670 </td>
   <td style="text-align:center;"> 27.40843 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita maxima </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> 123.9720248 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> 0.0020601 </td>
   <td style="text-align:center;"> 10.34886 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Cucurbita pepo </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> -26.7803739 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0096076 </td>
   <td style="text-align:center;"> 10.78700 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Dioscorea alata </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> 5.6181083 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> 0.0059156 </td>
   <td style="text-align:center;"> 23.84612 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Draba nivalis </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 34.11541 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Elaeis guineensis </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> 58.4634143 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> -0.0004549 </td>
   <td style="text-align:center;"> 46.94000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Eucalyptus grandis </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -52.1347009 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0159960 </td>
   <td style="text-align:center;"> 55.68998 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Glycine max </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -153.0590573 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0691302 </td>
   <td style="text-align:center;"> 48.12719 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium hirsutum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -92.0853380 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0647436 </td>
   <td style="text-align:center;"> 85.92276 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Gossypium raimondii </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -894.6366733 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0720796 </td>
   <td style="text-align:center;"> 57.63293 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Helianthus annuus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -92.8265579 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0291271 </td>
   <td style="text-align:center;"> 169.70400 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Hordeum vulgare </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -462.9124143 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.1904018 </td>
   <td style="text-align:center;"> 654.85949 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Juglans regia </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 47.43156 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus albus </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -21.5359754 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0291837 </td>
   <td style="text-align:center;"> 17.35603 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Lupinus angustifolius </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 27.89545 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Malus domestica </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -162.8282246 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0206064 </td>
   <td style="text-align:center;"> 38.38249 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Mangifera indica </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 17.87202 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Manihot esculenta </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> 486.9496847 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> -0.0099806 </td>
   <td style="text-align:center;"> 28.80229 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Momordica charantia </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 26.50729 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Nelumbo nucifera </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 97.84551 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza nivara </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -17.1467421 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0138083 </td>
   <td style="text-align:center;"> 28.16253 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Oryza sativa </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -132.8990325 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0236937 </td>
   <td style="text-align:center;"> 31.10379 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Panicum hallii </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -1922.1715758 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0865313 </td>
   <td style="text-align:center;"> 56.37928 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Phaseolus vulgaris </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -189.3509884 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0634829 </td>
   <td style="text-align:center;"> 46.80187 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus mume </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> 0.7466590 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> 0.0037906 </td>
   <td style="text-align:center;"> 24.98833 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Prunus persica </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -62.8750625 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0135736 </td>
   <td style="text-align:center;"> 28.21185 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Quercus sp </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 59.72765 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Raphanus sativus </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 38.78895 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sesamum indicum </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> exception </td>
   <td style="text-align:center;"> 12.6191473 </td>
   <td style="text-align:center;"> less </td>
   <td style="text-align:center;"> 0.0113517 </td>
   <td style="text-align:center;"> 16.40985 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Setaria italica </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -115.0067495 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0376572 </td>
   <td style="text-align:center;"> 44.58849 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum lycopersicum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -289.8449953 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.1352258 </td>
   <td style="text-align:center;"> 67.26872 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Solanum tuberosum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -101.5698020 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0529184 </td>
   <td style="text-align:center;"> 60.41812 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Sorghum bicolor </td>
   <td style="text-align:center;"> sub-telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -113.0771627 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.1995367 </td>
   <td style="text-align:center;"> 68.36450 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Theobroma cacao </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -189.3678093 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0415695 </td>
   <td style="text-align:center;"> 33.04562 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum aestivum </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -508.7577718 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0794583 </td>
   <td style="text-align:center;"> 667.82536 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum dicoccoides </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 582.40116 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Triticum urartu </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> NA </td>
   <td style="text-align:center;"> 665.94253 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vigna unguiculata </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -354.5383861 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0511767 </td>
   <td style="text-align:center;"> 43.04148 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Vitis vinifera </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> peak </td>
   <td style="text-align:center;"> -82.5057959 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.0211167 </td>
   <td style="text-align:center;"> 22.43032 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> Zea mays </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> telomere </td>
   <td style="text-align:center;"> -372.1453979 </td>
   <td style="text-align:center;"> more </td>
   <td style="text-align:center;"> 0.1188828 </td>
   <td style="text-align:center;"> 210.63381 </td>
  </tr>
</tbody>
</table>


```r
df = pattern_convergence
colnames(df) = c("Species", "Recombination pattern", "Gene density pattern", "RMSE difference", 
    "Homogeneization effect of gene distance", "Genetic shuffling difference", "Mean chromosome size")
write.xlsx(x = df, file = paste(wd, "tables/article_one_supp/tables_supplementary.xls", 
    sep = ""), sheetName = "S11_ConvergencePatterns", row.names = FALSE, append = TRUE)
rm(df)
```



```r
ggplot(data = pattern_convergence, aes(x = pattern, y = diffgeneticshuffl)) + geom_violin(aes(colour = pattern, 
    fill = pattern)) + geom_dotplot(aes(), binaxis = "y", stackdir = "center", dotsize = 0.5, 
    fill = 1) + scale_color_manual(values = c("black", colorpalette[2], colorpalette[1])) + 
    scale_fill_manual(values = alpha(c("black", colorpalette[2], colorpalette[1]), 
        0.5)) + xlab("Crossover pattern") + ylab("Difference in genetic shuffling") + 
    theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), 
        plot.title = element_text(color = "black", size = 14, face = "bold.italic", 
            hjust = 0.5), axis.title.x = element_text(color = "black", size = 14), 
        axis.title.y = element_text(color = "black", size = 14), axis.text = element_text(size = 14, 
            colour = "black"), legend.key = element_rect(fill = "white", size = 1), 
        legend.key.height = unit(1, "line"), legend.key.width = unit(2, "line"), 
        legend.text = element_text(size = 7, face = "italic"), legend.title = element_blank(), 
        legend.position = "none")
```

```
## Warning: Removed 17 rows containing non-finite values (stat_ydensity).
```

```
## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 17 rows containing non-finite values (stat_bindot).
```

![](/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-199-1.png)<!-- -->

```r
# plot(diffgeneticshuffl ~ log(meanchrsize), data = pattern_convergence)
```





```
## 
## Call:
## lm(formula = diff ~ chrsize, data = df)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.07294 -0.02894 -0.01314  0.01741  0.15455 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 3.274e-02  8.936e-03   3.663 0.000756 ***
## chrsize     1.792e-04  5.424e-05   3.303 0.002089 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.04847 on 38 degrees of freedom
## Multiple R-squared:  0.2231,	Adjusted R-squared:  0.2026 
## F-statistic: 10.91 on 1 and 38 DF,  p-value: 0.002089
```

<div class="figure" style="text-align: center">
<img src="/Users/tbrazier/Academic/PhD/diversity-determinants-recombination-landscapes-flowering-plants/reports/output/Results_final_files/figure-html/unnamed-chunk-200-1.png" alt="(ref:differences-shufflingrates)"  />
<p class="caption">(\#fig:unnamed-chunk-200)(ref:differences-shufflingrates)</p>
</div>

(ref:dist-shufflingrates) (ref:dist-shufflingrates)

#### Figure 11





# Packages & Session Info


```r
sessionInfo()
```

```
## R version 4.0.4 (2021-02-15)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] egg_0.4.5           ggrepel_0.9.1       ggtext_0.1.1       
##  [4] magick_2.7.2        gridExtra_2.3       xlsx_0.6.5         
##  [7] data.table_1.14.0   ggnewscale_0.4.5    treeio_1.14.3      
## [10] ggtree_2.4.1        MuMIn_1.43.17       pbmcapply_1.5.0    
## [13] reldist_1.6-6       lmerTest_3.1-3      plyr_1.8.6         
## [16] ape_5.4-1           car_3.0-10          carData_3.0-4      
## [19] kableExtra_1.3.4    phyr_1.1.0          INLA_21.02.23      
## [22] sjmisc_2.8.6        sjPlot_2.8.7        lme4_1.1-26        
## [25] Matrix_1.3-2        cowplot_1.1.1       pals_1.6           
## [28] ggplotify_0.0.5     RColorBrewer_1.1-2  reshape2_1.4.4     
## [31] ggpubr_0.4.0        purrr_0.3.4         rBLAST_0.99.2      
## [34] Biostrings_2.58.0   XVector_0.30.0      IRanges_2.24.1     
## [37] S4Vectors_0.28.1    BiocGenerics_0.36.0 rentrez_1.2.3      
## [40] dplyr_1.0.5         taxizedb_0.3.0      taxize_0.9.99      
## [43] gdata_2.18.0        stringr_1.4.0       MareyMap_1.3.6     
## [46] ggplot2_3.3.3       ade4_1.7-16         rstudioapi_0.13    
## [49] bookdown_0.21       knitr_1.38         
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.1          tidyselect_1.1.0    htmlwidgets_1.5.3  
##   [4] RSQLite_2.2.4       grid_4.0.4          hoardr_0.5.2       
##   [7] munsell_0.5.0       codetools_0.2-18    effectsize_0.4.4   
##  [10] statmod_1.4.35      withr_2.4.1         colorspace_2.0-0   
##  [13] highr_0.8           uuid_0.1-4          ggsignif_0.6.1     
##  [16] rJava_0.9-13        labeling_0.4.2      emmeans_1.5.4      
##  [19] conditionz_0.1.0    farver_2.1.0        bit64_4.0.5        
##  [22] coda_0.19-4         vctrs_0.3.8         generics_0.1.0     
##  [25] xfun_0.30           markdown_1.1        R6_2.5.0           
##  [28] cachem_1.0.4        reshape_0.8.8       gridGraphics_0.5-1 
##  [31] assertthat_0.2.1    scales_1.1.1        nnet_7.3-15        
##  [34] gtable_0.3.0        bold_1.1.0          rlang_1.0.2        
##  [37] systemfonts_1.0.1   splines_4.0.4       lazyeval_0.2.2     
##  [40] rstatix_0.7.0       dichromat_2.0-0     checkmate_2.0.0    
##  [43] broom_0.7.12        BiocManager_1.30.10 yaml_2.2.1         
##  [46] abind_1.4-5         modelr_0.1.8        backports_1.2.1    
##  [49] Hmisc_4.5-0         gridtext_0.1.4      tools_4.0.4        
##  [52] tcltk_4.0.4         ellipsis_0.3.2      jquerylib_0.1.3    
##  [55] Rcpp_1.0.6          base64enc_0.1-3     progress_1.2.2     
##  [58] zlibbioc_1.36.0     prettyunits_1.1.1   rpart_4.1-15       
##  [61] zoo_1.8-9           haven_2.3.1         cluster_2.1.1      
##  [64] crul_1.1.0          magrittr_2.0.1      openxlsx_4.2.3     
##  [67] mvtnorm_1.1-1       xlsxjars_0.6.1      patchwork_1.1.1    
##  [70] hms_1.0.0           evaluate_0.15       xtable_1.8-4       
##  [73] XML_3.99-0.5        rio_0.5.26          sjstats_0.18.1     
##  [76] jpeg_0.1-8.1        readxl_1.3.1        ggeffects_1.0.1    
##  [79] compiler_4.0.4      tibble_3.1.0        maps_3.3.0         
##  [82] crayon_1.4.1        minqa_1.2.4         htmltools_0.5.1.1  
##  [85] mgcv_1.8-34         Formula_1.2-4       aplot_0.0.6        
##  [88] tidyr_1.1.3         DBI_1.1.1           formatR_1.8        
##  [91] sjlabelled_1.1.7    dbplyr_2.1.1        MASS_7.3-53.1      
##  [94] rappdirs_0.3.3      boot_1.3-27         cli_3.2.0          
##  [97] insight_0.13.1      forcats_0.5.1       pkgconfig_2.0.3    
## [100] rvcheck_0.1.8       numDeriv_2016.8-1.1 foreign_0.8-81     
## [103] sp_1.4-5            xml2_1.3.3          foreach_1.5.1      
## [106] svglite_2.0.0       bslib_0.2.4         webshot_0.5.2      
## [109] estimability_1.3    rvest_1.0.0         digest_0.6.27      
## [112] parameters_0.12.0   httpcode_0.3.0      rmarkdown_2.7      
## [115] cellranger_1.1.0    tidytree_0.3.3      htmlTable_2.1.0    
## [118] curl_4.3            gtools_3.8.2        nloptr_1.2.2.2     
## [121] lifecycle_1.0.0     nlme_3.1-152        jsonlite_1.7.2     
## [124] mapproj_1.2.7       viridisLite_0.3.0   fansi_0.4.2        
## [127] pillar_1.7.0        lattice_0.20-41     fastmap_1.1.0      
## [130] httr_1.4.2          survival_3.2-9      glue_1.6.2         
## [133] bayestestR_0.8.2    zip_2.1.1           png_0.1-7          
## [136] iterators_1.0.13    bit_4.0.4           stringi_1.7.6      
## [139] sass_0.3.1          performance_0.7.0   blob_1.2.1         
## [142] latticeExtra_0.6-29 memoise_2.0.0
```

```r
# Citations for libraries
lapply(packages, function(x) citation(x))
```

```
## [[1]]
## 
## To cite package 'rstudioapi' in publications use:
## 
##   Kevin Ushey, JJ Allaire, Hadley Wickham and Gary Ritchie (2020).
##   rstudioapi: Safely Access the RStudio API. R package version 0.13.
##   https://CRAN.R-project.org/package=rstudioapi
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {rstudioapi: Safely Access the RStudio API},
##     author = {Kevin Ushey and JJ Allaire and Hadley Wickham and Gary Ritchie},
##     year = {2020},
##     note = {R package version 0.13},
##     url = {https://CRAN.R-project.org/package=rstudioapi},
##   }
## 
## 
## [[2]]
## 
## To cite ade4 in publications use:
## 
## Dray S, Dufour A (2007). "The ade4 Package: Implementing the Duality
## Diagram for Ecologists." _Journal of Statistical Software_, *22*(4),
## 1-20. doi: 10.18637/jss.v022.i04 (URL:
## https://doi.org/10.18637/jss.v022.i04).
## 
## Bougeard S, Dray S (2018). "Supervised Multiblock Analysis in R with
## the ade4 Package." _Journal of Statistical Software_, *86*(1), 1-17.
## doi: 10.18637/jss.v086.i01 (URL:
## https://doi.org/10.18637/jss.v086.i01).
## 
## Chessel D, Dufour A, Thioulouse J (2004). "The ade4 Package - I:
## One-Table Methods." _R News_, *4*(1), 5-10. <URL:
## https://cran.r-project.org/doc/Rnews/>.
## 
## Dray S, Dufour A, Chessel D (2007). "The ade4 Package - II: Two-Table
## and K-Table Methods." _R News_, *7*(2), 47-52. <URL:
## https://cran.r-project.org/doc/Rnews/>.
## 
## Thioulouse J, Dray S, Dufour A, Siberchicot A, Jombart T, Pavoine S
## (2018). _Multivariate Analysis of Ecological Data with ade4_. Springer.
## doi: 10.1007/978-1-4939-8850-1 (URL:
## https://doi.org/10.1007/978-1-4939-8850-1).
## 
## To see these entries in BibTeX format, use 'print(<citation>,
## bibtex=TRUE)', 'toBibtex(.)', or set
## 'options(citation.bibtex.max=999)'.
## 
## 
## [[3]]
## 
## To cite ggplot2 in publications, please use:
## 
##   H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
##   Springer-Verlag New York, 2016.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Book{,
##     author = {Hadley Wickham},
##     title = {ggplot2: Elegant Graphics for Data Analysis},
##     publisher = {Springer-Verlag New York},
##     year = {2016},
##     isbn = {978-3-319-24277-4},
##     url = {https://ggplot2.tidyverse.org},
##   }
## 
## 
## [[4]]
## 
## To cite package 'MareyMap' in publications use:
## 
##   Aurélie Siberchicot, Clément Rezvoy, Delphine Charif, Laurent Gueguen
##   and Gabriel Marais (2020). MareyMap: Estimation of Meiotic
##   Recombination Rates Using Marey Maps. R package version 1.3.6.
##   https://CRAN.R-project.org/package=MareyMap
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {MareyMap: Estimation of Meiotic Recombination Rates Using Marey Maps},
##     author = {Aurélie Siberchicot and Clément Rezvoy and Delphine Charif and Laurent Gueguen and Gabriel Marais},
##     year = {2020},
##     note = {R package version 1.3.6},
##     url = {https://CRAN.R-project.org/package=MareyMap},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
## 
## 
## [[5]]
## 
## To cite package 'stringr' in publications use:
## 
##   Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for
##   Common String Operations. R package version 1.4.0.
##   https://CRAN.R-project.org/package=stringr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {stringr: Simple, Consistent Wrappers for Common String Operations},
##     author = {Hadley Wickham},
##     year = {2019},
##     note = {R package version 1.4.0},
##     url = {https://CRAN.R-project.org/package=stringr},
##   }
## 
## 
## [[6]]
## 
## To cite package 'gdata' in publications use:
## 
##   Gregory R. Warnes, Ben Bolker, Gregor Gorjanc, Gabor Grothendieck,
##   Ales Korosec, Thomas Lumley, Don MacQueen, Arni Magnusson, Jim Rogers
##   and others (2017). gdata: Various R Programming Tools for Data
##   Manipulation. R package version 2.18.0.
##   https://CRAN.R-project.org/package=gdata
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {gdata: Various R Programming Tools for Data Manipulation},
##     author = {Gregory R. Warnes and Ben Bolker and Gregor Gorjanc and Gabor Grothendieck and Ales Korosec and Thomas Lumley and Don MacQueen and Arni Magnusson and Jim Rogers and others},
##     year = {2017},
##     note = {R package version 2.18.0},
##     url = {https://CRAN.R-project.org/package=gdata},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
## 
## 
## [[7]]
## 
## To cite taxize in publications use at least the first, if not both:
## 
##   Scott Chamberlain and Eduard Szocs (2013). taxize - taxonomic search
##   and retrieval in R. F1000Research, 2:191. URL:
##   https://f1000research.com/articles/2-191/v2
## 
##   Scott Chamberlain, Eduard Szoecs, Zachary Foster, Zebulun Arendsee,
##   Carl Boettiger, Karthik Ram, Ignasi Bartomeus, John Baumgartner,
##   James O'Donnell, Jari Oksanen, Bastian Greshake Tzovaras, Philippe
##   Marchand, Vinh Tran, Maëlle Salmon, Gaopeng Li, and Matthias Grenié.
##   (2020) taxize: Taxonomic information from around the web. R package
##   version 0.9.98. https://github.com/ropensci/taxize
## 
## To see these entries in BibTeX format, use 'print(<citation>,
## bibtex=TRUE)', 'toBibtex(.)', or set
## 'options(citation.bibtex.max=999)'.
## 
## 
## [[8]]
## 
## To cite package 'taxizedb' in publications use:
## 
##   Scott Chamberlain and Zebulun Arendsee (2021). taxizedb: Tools for
##   Working with 'Taxonomic' Databases. R package version 0.3.0.
##   https://CRAN.R-project.org/package=taxizedb
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {taxizedb: Tools for Working with 'Taxonomic' Databases},
##     author = {Scott Chamberlain and Zebulun Arendsee},
##     year = {2021},
##     note = {R package version 0.3.0},
##     url = {https://CRAN.R-project.org/package=taxizedb},
##   }
## 
## 
## [[9]]
## 
## To cite package 'dplyr' in publications use:
## 
##   Hadley Wickham, Romain François, Lionel Henry and Kirill Müller
##   (2021). dplyr: A Grammar of Data Manipulation. R package version
##   1.0.5. https://CRAN.R-project.org/package=dplyr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {dplyr: A Grammar of Data Manipulation},
##     author = {Hadley Wickham and Romain François and Lionel Henry and Kirill Müller},
##     year = {2021},
##     note = {R package version 1.0.5},
##     url = {https://CRAN.R-project.org/package=dplyr},
##   }
## 
## 
## [[10]]
## 
## To cite rentrez in publications use:
## 
##   Winter, D. J. (2017) rentrez: an R package for the NCBI eUtils API
##   The R Journal 9(2):520-526
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {{rentrez}: an R package for the NCBI eUtils API},
##     author = {David J. Winter},
##     journal = {The R Journal},
##     year = {2017},
##     volume = {9},
##     issue = {2},
##     pages = {520--526},
##   }
## 
## 
## [[11]]
## 
## To cite package 'rBLAST' in publications use:
## 
##   Michael Hahsler and Anurag Nagar (2019). rBLAST: Interface to the
##   Basic Local Alignment Search Tool (BLAST). R package version 0.99.2.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {rBLAST: Interface to the Basic Local Alignment Search Tool (BLAST)},
##     author = {Michael Hahsler and Anurag Nagar},
##     year = {2019},
##     note = {R package version 0.99.2},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
## 
## 
## [[12]]
## 
## To cite package 'Biostrings' in publications use:
## 
##   H. Pagès, P. Aboyoun, R. Gentleman and S. DebRoy (2021). Biostrings:
##   Efficient manipulation of biological strings. R package version
##   2.58.0. https://bioconductor.org/packages/Biostrings
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {Biostrings: Efficient manipulation of biological strings},
##     author = {H. Pagès and P. Aboyoun and R. Gentleman and S. DebRoy},
##     year = {2021},
##     note = {R package version 2.58.0},
##     url = {https://bioconductor.org/packages/Biostrings},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
## 
## 
## [[13]]
## 
## To cite package 'purrr' in publications use:
## 
##   Lionel Henry and Hadley Wickham (2020). purrr: Functional Programming
##   Tools. R package version 0.3.4.
##   https://CRAN.R-project.org/package=purrr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {purrr: Functional Programming Tools},
##     author = {Lionel Henry and Hadley Wickham},
##     year = {2020},
##     note = {R package version 0.3.4},
##     url = {https://CRAN.R-project.org/package=purrr},
##   }
## 
## 
## [[14]]
## 
## To cite package 'ggpubr' in publications use:
## 
##   Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication
##   Ready Plots. R package version 0.4.0.
##   https://CRAN.R-project.org/package=ggpubr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {ggpubr: 'ggplot2' Based Publication Ready Plots},
##     author = {Alboukadel Kassambara},
##     year = {2020},
##     note = {R package version 0.4.0},
##     url = {https://CRAN.R-project.org/package=ggpubr},
##   }
## 
## 
## [[15]]
## 
## To cite reshape2 in publications use:
## 
##   Hadley Wickham (2007). Reshaping Data with the reshape Package.
##   Journal of Statistical Software, 21(12), 1-20. URL
##   http://www.jstatsoft.org/v21/i12/.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Reshaping Data with the {reshape} Package},
##     author = {Hadley Wickham},
##     journal = {Journal of Statistical Software},
##     year = {2007},
##     volume = {21},
##     number = {12},
##     pages = {1--20},
##     url = {http://www.jstatsoft.org/v21/i12/},
##   }
## 
## 
## [[16]]
## 
## To cite package 'RColorBrewer' in publications use:
## 
##   Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package
##   version 1.1-2. https://CRAN.R-project.org/package=RColorBrewer
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {RColorBrewer: ColorBrewer Palettes},
##     author = {Erich Neuwirth},
##     year = {2014},
##     note = {R package version 1.1-2},
##     url = {https://CRAN.R-project.org/package=RColorBrewer},
##   }
## 
## 
## [[17]]
## 
## To cite package 'ggplotify' in publications use:
## 
##   Guangchuang Yu (2020). ggplotify: Convert Plot to 'grob' or 'ggplot'
##   Object. R package version 0.0.5.
##   https://CRAN.R-project.org/package=ggplotify
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {ggplotify: Convert Plot to 'grob' or 'ggplot' Object},
##     author = {Guangchuang Yu},
##     year = {2020},
##     note = {R package version 0.0.5},
##     url = {https://CRAN.R-project.org/package=ggplotify},
##   }
## 
## 
## [[18]]
## 
## To cite package 'pals' in publications use:
## 
##   Kevin Wright (2019). pals: Color Palettes, Colormaps, and Tools to
##   Evaluate Them. R package version 1.6.
##   https://CRAN.R-project.org/package=pals
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {pals: Color Palettes, Colormaps, and Tools to Evaluate Them},
##     author = {Kevin Wright},
##     year = {2019},
##     note = {R package version 1.6},
##     url = {https://CRAN.R-project.org/package=pals},
##   }
## 
## 
## [[19]]
## 
## To cite package 'cowplot' in publications use:
## 
##   Claus O. Wilke (2020). cowplot: Streamlined Plot Theme and Plot
##   Annotations for 'ggplot2'. R package version 1.1.1.
##   https://CRAN.R-project.org/package=cowplot
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {cowplot: Streamlined Plot Theme and Plot Annotations for 'ggplot2'},
##     author = {Claus O. Wilke},
##     year = {2020},
##     note = {R package version 1.1.1},
##     url = {https://CRAN.R-project.org/package=cowplot},
##   }
## 
## 
## [[20]]
## 
## To cite lme4 in publications use:
## 
##   Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015).
##   Fitting Linear Mixed-Effects Models Using lme4. Journal of
##   Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Fitting Linear Mixed-Effects Models Using {lme4}},
##     author = {Douglas Bates and Martin M{\"a}chler and Ben Bolker and Steve Walker},
##     journal = {Journal of Statistical Software},
##     year = {2015},
##     volume = {67},
##     number = {1},
##     pages = {1--48},
##     doi = {10.18637/jss.v067.i01},
##   }
## 
## 
## [[21]]
## 
## Lüdecke D (2021). _sjPlot: Data Visualization for Statistics in Social
## Science_. R package version 2.8.7, <URL:
## https://CRAN.R-project.org/package=sjPlot>.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {sjPlot: Data Visualization for Statistics in Social Science},
##     author = {Daniel Lüdecke},
##     year = {2021},
##     note = {R package version 2.8.7},
##     url = {https://CRAN.R-project.org/package=sjPlot},
##   }
## 
## 
## [[22]]
## 
## Lüdecke D (2018). "sjmisc: Data and Variable Transformation Functions."
## _Journal of Open Source Software_, *3*(26), 754. doi:
## 10.21105/joss.00754 (URL: https://doi.org/10.21105/joss.00754).
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {sjmisc: Data and Variable Transformation Functions.},
##     volume = {3},
##     doi = {10.21105/joss.00754},
##     number = {26},
##     journal = {Journal of Open Source Software},
##     author = {Daniel Lüdecke},
##     year = {2018},
##     pages = {754},
##   }
## 
## 
## [[23]]
## 
## To cite package 'phyr' in publications use:
## 
##   Anthony Ives, Russell Dinnage, Lucas A. Nell, Matthew Helmus and
##   Daijiang Li (2020). phyr: Model Based Phylogenetic Analysis. R
##   package version 1.1.0. https://CRAN.R-project.org/package=phyr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {phyr: Model Based Phylogenetic Analysis},
##     author = {Anthony Ives and Russell Dinnage and Lucas A. Nell and Matthew Helmus and Daijiang Li},
##     year = {2020},
##     note = {R package version 1.1.0},
##     url = {https://CRAN.R-project.org/package=phyr},
##   }
## 
## 
## [[24]]
## 
## To cite package 'kableExtra' in publications use:
## 
##   Hao Zhu (2021). kableExtra: Construct Complex Table with 'kable' and
##   Pipe Syntax. R package version 1.3.4.
##   https://CRAN.R-project.org/package=kableExtra
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {kableExtra: Construct Complex Table with 'kable' and Pipe Syntax},
##     author = {Hao Zhu},
##     year = {2021},
##     note = {R package version 1.3.4},
##     url = {https://CRAN.R-project.org/package=kableExtra},
##   }
## 
## 
## [[25]]
## 
## To cite the car package in publications use:
## 
##   John Fox and Sanford Weisberg (2019). An {R} Companion to Applied
##   Regression, Third Edition. Thousand Oaks CA: Sage. URL:
##   https://socialsciences.mcmaster.ca/jfox/Books/Companion/
## 
## A BibTeX entry for LaTeX users is
## 
##   @Book{,
##     title = {An {R} Companion to Applied Regression},
##     edition = {Third},
##     author = {John Fox and Sanford Weisberg},
##     year = {2019},
##     publisher = {Sage},
##     address = {Thousand Oaks {CA}},
##     url = {https://socialsciences.mcmaster.ca/jfox/Books/Companion/},
##   }
## 
## 
## [[26]]
## 
## To cite ape in a publication please use:
## 
##   Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern
##   phylogenetics and evolutionary analyses in R. Bioinformatics 35:
##   526-528.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {ape 5.0: an environment for modern phylogenetics and evolutionary analyses in {R}},
##     author = {E. Paradis and K. Schliep},
##     journal = {Bioinformatics},
##     year = {2019},
##     volume = {35},
##     pages = {526-528},
##   }
## 
## As ape is evolving quickly, you may want to cite also its version
## number (found with 'library(help = ape)' or 'packageVersion("ape")').
## 
## 
## [[27]]
## 
## To cite plyr in publications use:
## 
##   Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data
##   Analysis. Journal of Statistical Software, 40(1), 1-29. URL
##   http://www.jstatsoft.org/v40/i01/.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {The Split-Apply-Combine Strategy for Data Analysis},
##     author = {Hadley Wickham},
##     journal = {Journal of Statistical Software},
##     year = {2011},
##     volume = {40},
##     number = {1},
##     pages = {1--29},
##     url = {http://www.jstatsoft.org/v40/i01/},
##   }
## 
## 
## [[28]]
## 
## To cite lmerTest in publications use:
## 
## Kuznetsova A, Brockhoff PB, Christensen RHB (2017). "lmerTest Package:
## Tests in Linear Mixed Effects Models." _Journal of Statistical
## Software_, *82*(13), 1-26. doi: 10.18637/jss.v082.i13 (URL:
## https://doi.org/10.18637/jss.v082.i13).
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {{lmerTest} Package: Tests in Linear Mixed Effects Models},
##     author = {Alexandra Kuznetsova and Per B. Brockhoff and Rune H. B. Christensen},
##     journal = {Journal of Statistical Software},
##     year = {2017},
##     volume = {82},
##     number = {13},
##     pages = {1--26},
##     doi = {10.18637/jss.v082.i13},
##   }
## 
## 
## [[29]]
## 
## If you are using ideas from the 'Relative Distribution Methods' book
## for research that will be published, we request that you acknowledge
## this by citing the book as shown below.
## 
## If you use this 'reldist' package, please also use the second citation
## below.
## 
## For BibTeX format, use toBibtex(citation("reldist")).
## 
##   Mark S. Handcock and Martina Morris (1999) Relative Distribution
##   Methods in the Social Sciences. Springer, New York, ISBN
##   0-387-98778-9. URL http://www.stat.ucla.edu/~handcock/RelDist
## 
##   Mark S. Handcock (2016), Relative Distribution Methods. Version
##   1.6-6. Project home page at
##   http://www.stat.ucla.edu/~handcock/RelDist. URL
##   https://CRAN.R-project.org/package=reldist.
## 
## We have invested a lot of time and effort in creating the 'reldist'
## package for use by other researchers. Please cite it in all papers
## where it is used.
## To see these entries in BibTeX format, use 'print(<citation>,
## bibtex=TRUE)', 'toBibtex(.)', or set
## 'options(citation.bibtex.max=999)'.
## 
## 
## [[30]]
## 
## To cite package 'pbmcapply' in publications use:
## 
##   Kevin Kuang, Quyu Kong and Francesco Napolitano (2019). pbmcapply:
##   Tracking the Progress of Mc*pply with Progress Bar. R package version
##   1.5.0. https://CRAN.R-project.org/package=pbmcapply
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {pbmcapply: Tracking the Progress of Mc*pply with Progress Bar},
##     author = {Kevin Kuang and Quyu Kong and Francesco Napolitano},
##     year = {2019},
##     note = {R package version 1.5.0},
##     url = {https://CRAN.R-project.org/package=pbmcapply},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
## 
## 
## [[31]]
## 
## To cite package 'MuMIn' in publications use:
## 
##   Kamil Bartoń (2020). MuMIn: Multi-Model Inference. R package version
##   1.43.17. https://CRAN.R-project.org/package=MuMIn
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {MuMIn: Multi-Model Inference},
##     author = {Kamil Bartoń},
##     year = {2020},
##     note = {R package version 1.43.17},
##     url = {https://CRAN.R-project.org/package=MuMIn},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
## 
## 
## [[32]]
## 
## To cite ggtree in publications use:
## 
##   Guangchuang Yu. Using ggtree to visualize data on tree-like
##   structures. Current Protocols in Bioinformatics, 2020, 69:e96. doi:
##   10.1002/cpbi.96
## 
##   Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods
##   for mapping and visualizing associated data on phylogeny using
##   ggtree. Molecular Biology and Evolution 2018, 35(2):3041-3043. doi:
##   10.1093/molbev/msy194
## 
##   Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk
##   Lam. ggtree: an R package for visualization and annotation of
##   phylogenetic trees with their covariates and other associated data.
##   Methods in Ecology and Evolution 2017, 8(1):28-36.
##   doi:10.1111/2041-210X.12628
## 
## To see these entries in BibTeX format, use 'print(<citation>,
## bibtex=TRUE)', 'toBibtex(.)', or set
## 'options(citation.bibtex.max=999)'.
## 
## 
## [[33]]
## 
## To cite treeio in publications use:
## 
##   LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
##   Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
##   for phylogenetic tree input and output with richly annotated and
##   associated data. Molecular Biology and Evolution 2020, 37(2):599-603.
##   doi: 10.1093/molbev/msz240
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {treeio: an R package for phylogenetic tree input and output with richly annotated and associated data.},
##     author = {Li-Gen Wang and Tommy Tsan-Yuk Lam and Shuangbin Xu and Zehan Dai and Lang Zhou and Tingze Feng and Pingfan Guo and Casey W. Dunn and Bradley R. Jones and Tyler Bradley and Huachen Zhu and Yi Guan and Yong Jiang and Guangchuang Yu},
##     year = {2020},
##     journal = {Molecular Biology and Evolution},
##     volume = {37},
##     issue = {2},
##     pages = {599-603},
##     doi = {10.1093/molbev/msz240},
##   }
## 
## 
## [[34]]
## 
## To cite package 'ggnewscale' in publications use:
## 
##   Elio Campitelli (2021). ggnewscale: Multiple Fill and Colour Scales
##   in 'ggplot2'. R package version 0.4.5.
##   https://CRAN.R-project.org/package=ggnewscale
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {ggnewscale: Multiple Fill and Colour Scales in 'ggplot2'},
##     author = {Elio Campitelli},
##     year = {2021},
##     note = {R package version 0.4.5},
##     url = {https://CRAN.R-project.org/package=ggnewscale},
##   }
## 
## 
## [[35]]
## 
## To cite package 'dplyr' in publications use:
## 
##   Hadley Wickham, Romain François, Lionel Henry and Kirill Müller
##   (2021). dplyr: A Grammar of Data Manipulation. R package version
##   1.0.5. https://CRAN.R-project.org/package=dplyr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {dplyr: A Grammar of Data Manipulation},
##     author = {Hadley Wickham and Romain François and Lionel Henry and Kirill Müller},
##     year = {2021},
##     note = {R package version 1.0.5},
##     url = {https://CRAN.R-project.org/package=dplyr},
##   }
## 
## 
## [[36]]
## 
## To cite package 'data.table' in publications use:
## 
##   Matt Dowle and Arun Srinivasan (2021). data.table: Extension of
##   `data.frame`. R package version 1.14.0.
##   https://CRAN.R-project.org/package=data.table
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {data.table: Extension of `data.frame`},
##     author = {Matt Dowle and Arun Srinivasan},
##     year = {2021},
##     note = {R package version 1.14.0},
##     url = {https://CRAN.R-project.org/package=data.table},
##   }
## 
## 
## [[37]]
## 
## To cite package 'xlsx' in publications use:
## 
##   Adrian Dragulescu and Cole Arendt (2020). xlsx: Read, Write, Format
##   Excel 2007 and Excel 97/2000/XP/2003 Files. R package version 0.6.5.
##   https://CRAN.R-project.org/package=xlsx
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {xlsx: Read, Write, Format Excel 2007 and Excel 97/2000/XP/2003 Files},
##     author = {Adrian Dragulescu and Cole Arendt},
##     year = {2020},
##     note = {R package version 0.6.5},
##     url = {https://CRAN.R-project.org/package=xlsx},
##   }
## 
## 
## [[38]]
## 
## To cite package 'gridExtra' in publications use:
## 
##   Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid"
##   Graphics. R package version 2.3.
##   https://CRAN.R-project.org/package=gridExtra
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {gridExtra: Miscellaneous Functions for "Grid" Graphics},
##     author = {Baptiste Auguie},
##     year = {2017},
##     note = {R package version 2.3},
##     url = {https://CRAN.R-project.org/package=gridExtra},
##   }
## 
## 
## [[39]]
## 
## To cite package 'magick' in publications use:
## 
##   Jeroen Ooms (2021). magick: Advanced Graphics and Image-Processing in
##   R. R package version 2.7.2. https://CRAN.R-project.org/package=magick
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {magick: Advanced Graphics and Image-Processing in R},
##     author = {Jeroen Ooms},
##     year = {2021},
##     note = {R package version 2.7.2},
##     url = {https://CRAN.R-project.org/package=magick},
##   }
## 
## 
## [[40]]
## 
## To cite package 'ggtext' in publications use:
## 
##   Claus O. Wilke (2020). ggtext: Improved Text Rendering Support for
##   'ggplot2'. R package version 0.1.1.
##   https://CRAN.R-project.org/package=ggtext
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {ggtext: Improved Text Rendering Support for 'ggplot2'},
##     author = {Claus O. Wilke},
##     year = {2020},
##     note = {R package version 0.1.1},
##     url = {https://CRAN.R-project.org/package=ggtext},
##   }
## 
## 
## [[41]]
## 
## To cite package 'ggrepel' in publications use:
## 
##   Kamil Slowikowski (2021). ggrepel: Automatically Position
##   Non-Overlapping Text Labels with 'ggplot2'. R package version 0.9.1.
##   https://CRAN.R-project.org/package=ggrepel
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {ggrepel: Automatically Position Non-Overlapping Text Labels with
## 'ggplot2'},
##     author = {Kamil Slowikowski},
##     year = {2021},
##     note = {R package version 0.9.1},
##     url = {https://CRAN.R-project.org/package=ggrepel},
##   }
## 
## 
## [[42]]
## 
## To cite package 'egg' in publications use:
## 
##   Baptiste Auguie (2019). egg: Extensions for 'ggplot2': Custom Geom,
##   Custom Themes, Plot Alignment, Labelled Panels, Symmetric Scales, and
##   Fixed Panel Size. R package version 0.4.5.
##   https://CRAN.R-project.org/package=egg
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {egg: Extensions for 'ggplot2': Custom Geom, Custom Themes, Plot
## Alignment, Labelled Panels, Symmetric Scales, and Fixed Panel
## Size},
##     author = {Baptiste Auguie},
##     year = {2019},
##     note = {R package version 0.4.5},
##     url = {https://CRAN.R-project.org/package=egg},
##   }
```

# References



