# Assignment 4_Bioinformatics Programming Challenges

To obtain the .txt file with the orthologue pairs between species A. thaliana and S. pombe, it is necessary to execute the following command:

**$ ruby assignment4.rb ./BLAST_databases/arabidopsis_th.fa ./BLAST_databases/s_pombe.fa**

## WARNING: since the computation time needed to perform the search for all ortholog pairs is quite high (in principle it should take a couple of hours, but in my case the execution of the code was always interrupted after several minutes, without any apparent reason), it is recommended to use the following command instead:

**$ ruby assignment4_test.rb ./BLAST_databases/arabidopsis_th.fa ./BLAST_databases/s_pombe.fa**

As some kind of problem has arisen that prevented the execution of the final script, 3 different .txt files have been generated, testing different number of queries (100, 1000 and 3000, respectively). 
This is not as optimal as if the complete result had been obtained, but *the greater the number of queries, the greater the similarity to a complete analysis*.

The principal BLAST parameters on which the search for orthologues has been based are the following:

* **E-value** --> threshold value stablished at **1×10−6**
* Detection of orthologues as best reciprocal hits with **Smith-Waterman final alignment (-F “m S”)**, giving both the highest number of orthologs and the minimal error rates.
* It has not been included, but it would have been interesting considering a **query coverage** of at least **50%** (as it is widely recommended in all references presented).

## Bibliograpy:
* https://doi.org/10.1093/bioinformatics/btm585
* https://dx.doi.org/10.1371%2Fjournal.pone.0101850
* https://doi.org/10.1186/s12864-020-07132-6

## BONUS (1%):

**Expanding the analysis of the putative orthologues discovered:**

* One interesting way would be **including a third species** (ideally related to Arabidopsis and S. pombe.), and perform some **clustering method** such as COG search (cluster of orthologous genes), based on *best reciprocal hits*.

* There are multiple *Tools and Methodologies* that are currently being developped in this field, such as **Pipelines for building Ortholog Data Sets**, **Comprehensive Visual Exploration in Orthology and Paralogy Analysis**, and many more examples, as *this article* shows (https://doi.org/10.3389/fgene.2017.00165)

