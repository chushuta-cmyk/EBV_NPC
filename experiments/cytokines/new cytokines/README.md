# New cytokine analysis

Run the analysis with:

```bash
Rscript analyze_new_cytokines.R
```

The workbook only identifies samples as 1 to 6. The script therefore plots them as separate conditions. To combine biological replicates, add `sample_map.csv` with these columns:

```text
Sample,Condition
1,DMSO
2,DMSO
3,DMSO
4,TPA
5,TPA
6,TPA
```

Change the example mapping to the true experimental design before using grouped means or error bars.
