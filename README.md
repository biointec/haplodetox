# HaploDetox
HaploDetox determines strain-specific de Bruijn graphs through accurate determination of strain-aware node and arc multiplicities in a mixed-sample de Bruijn graph using Conditional Random Fields.

A short paper describing HaploDetox has been submitted to CIBB 2023: 18th Conference on Computational Intelligence Methods for Bioinformatics & Biostatistics.

HaploDetox is an extension of [Detox](https://github.com/biointec/detox), our method for multiplicity estimation in de Bruijn graphs from isolate genomes.

## When to use HaploDetox

HaploDetox is a method to clean and split a de Bruijn graph into several strain-specific de Bruijn graphs.
It is aimed at **bacterial data**. 
We assume you are inputting a set of reads that have already been clustered by species, but are still a combination of a **small number of strains**.

Starting from a metagenomic dataset, different approaches exist to cluster reads by species. For example, you could map the reads to a reference database, use a k-mer-based clustering strategy, or map the reads to Metagenome Assembled Genomes (MAGs) after running a metagenomic assembler that is not strain-aware.

We assume the following about the data:
* The number of strains is small.

   When the number of strains is not passed to HaploDetox, we only test whether the number of strains present is 2, 3, or 4.  
   You can choose to determine the number of strains with an external method and pass this to HaploDetox, but number of strains $>4$ has not been benchmarked.

* The strains are very similar and a large part of their genome sequence is shared by all strains.

   We fit a model to the k-mer spectrum and this is an important model specification that is used during model training.
   More specifically, your k-mer spectrum should look like shown below (for two or three strains, other numbers of strains are possible). You will see several small peaks at the average coverage for individual strains (strain-specific multiplicity 1 in one strain and 0 in the others), and a large peak at a higher average coverage of the shared k-mers (strain-specific multiplicity 1 in all strains). The fitted mixture model has more components than these, but the unique-strain coverage and all-strain coverage peaks are usually the best visible. Note that the peaks corresponding to individual strains might be much smaller than shown here if strain divergence is very low, this should not be an issue. It is however important that the largest peak should correspond to the shared k-mers. We might extend HaploDetox to include other model assumptions in the future.

|         | Histogram | Fitted mixture model |
|---------|-----------|----------------------|
|2 strains| ![alt text](https://github.com/biointec/haplodetox/blob/main/example/spectra/2strains_histogram.png) | ![alt text](https://github.com/biointec/haplodetox/blob/main/example/spectra/2strains_modelfit.png) |
|3 strains| ![alt text](https://github.com/biointec/haplodetox/blob/main/example/spectra/3strains_spectrum.png) | ![alt text](https://github.com/biointec/haplodetox/blob/main/example/spectra/3strains_modelfit.png) |

## Prerequisites
To build HaploDetox you need the following software or library packages:

  * CMake version 2.6.3 or higher
  * GCC version 4.7 or higher
  * Google's [sparsehash](https://github.com/sparsehash/sparsehash)
  * [LAPACK](https://netlib.org/lapack/)
  * ZLIB (optional)
  * recent [Boost](https://www.boost.org/) libraries (for approximate inference)
  * [GMP library](https://gmplib.org/) (for approximate inference)

For approximate inference computations, HaploDetox relies on an adaptation of the Loopy Belief Propagation implementation of the [libDAI](https://bitbucket.org/jorism/libdai/src/master/) C++ library. 
The classes from libDAI that are necessary for Loopy Belief Propagation with our own additions of more recent message passing schemes are provided in this repository.

To use HaploDetox you also need to install and compile [BCALM 2](https://github.com/GATB/bcalm).

## Compilation

Clone the HaploDetox Github repository:

```bash
git clone https://github.com/biointec/haplodetox.git
```

Next, the C++ code can be compiled as follows:

```bash
cd haplodetox
mkdir build
cd build
cmake ..
make
```

After compilation, the `haplodetox` binary can be found in:

```bash
haplodetox/build/src/haplodetox
```

## Usage

We use BCALM 2 to build a compacted de Bruijn graph from sequencing data. HaploDetox determines strain-specific multiplicities for all nodes and arcs in this graph, cleans the graph, recomputes multiplicities, and finally uses the strain-specific multiplicities to phase the de Bruijn graph into several strain-specific de Bruijn graphs.
### 1. Generating the de Bruijn graph using BCALM 2
The first step is to prepare a manifest file, e.g. named `reads.mf`, in which the input FASTQ files are listed. For example, suppose the reads are stored in two FASTQ files named `reads_1.fastq` and `reads_2.fastq`, then the `reads.mf` file simply lists these input files, one file per line:

```
reads_1.fastq
reads_2.fastq
```
The set of reads must be in FASTQ (`.fq`/`.fastq`) format.

Using BCALM 2, one can then create a de Bruijn graph as follows:

```bash
bcalm -in reads.mf -kmer-size 21 -abundance-min 2
```

This process is usually very fast for bacterial genomes. The command line option `-abundance-min 2` instructs BCALM 2 to only retain k-mers that occur more than once in the input reads. BCALM 2's default value for `-abundance-min` is 2. We recommend using this value. However, in case of low sequencing depth you might consider using `-abundance-min 1` in BCALM and relying on HaploDetox' graph cleaning procedure alone to remove spurious k-mers. 
Setting this parameter to a larger value will remove a large fraction of the sequencing errors already, at the cost of possibly removing true k-mers.

When BCALM 2 is finished, it produces the file `reads.unitigs.fa`, in which the de Bruijn graph's unitig node sequences and their connectivity are encoded in a FASTA-format. This file serves as input for HaploDetox.

### 2. Retrieving phased contigs by inferring the strain aware node/arc multiplicities using HaploDetox

The pipeline takes as input a de Bruijn graph in the FASTA-format of BCALM 2 (`reads.unitigs.fa`) as well as the original read set (`reads.mf`). HaploDetox can then be run as follows:

```bash
haplodetox -mm-nstrains 2 -use-qual -abundance-min 2 reads.unitigs.fa reads.mf
```

The two mandatory input arguments are `reads.unitigs.fa` and `reads.mf`, in that order. All other flags or command line options should be listed prior to these two input arguments. 

If you have a prior estimate of the number of strains present in the data, you should pass this value with the `-mm-nstrains` parameter. If the parameter is not passed, HaploDetox trains a mixture model of the k-mer spectrum for 2, 3, and 4 strains and selects the number of strains for which the model fit is best.

The flag `-use-qual` instructs HaploDetox to weigh the k-mer counts using the Phred quality scores. When BCALM 2 was run using the `-abundance-min` option, it is best to provide this information also to HaploDetox. When HaploDetox is fitting the error model to the data, it assumes all k-mers that occur fewer times than the value provided by `abundance-min` are missing. It takes this into account when fitting the model to the data. If you use HaploDetox with unweighted k-mer counts, the value of `-abundance-min` can be automatically inferred.

Upon completion, the estimated node and arc multiplicities are stored in files `estmult.node` and `estmult.arc`, respectively, while the files `strainx.phased.fasta` contain the contigs of the de Bruijn graph that represents strain x in fasta format. There will be one file for each strain.
Note that HaploDetox does not perform a full haplotype-aware assembly. These contigs are much shorter than a typical assembly as no repeat resolution has been performed.

The file `estmult.node` contains the following fields (tab-separated):
* The node identifier (integer value from 1 to #nodes)
* The estimated multiplicity in strain 1 m1 (non-negative integer value)
* ...
* The estimated multiplicity in strain S mS (non-negative integer value)
* The log-probability _log[ P(mult. is (m1, ..., mS) ]_
* The average node coverage (weighed by the Phred quality use if the `-use-qual` flag was used)
* The length of the node expressed as the number of k-mers.
* The entropy of the assignment

The file `estmult.edge` contains the following fields (tab-separated):

* The node identifier of the node from which the arc is originating
* The node identifier of the destination node
* The estimated multiplicity m1 (non-negative integer value)
* ...
* The estimated multiplicity mS (non-negative integer value)
* The log-probability _log[ P(mult. is (m1, ..., mS)) ]_
* The edge coverage (weighed by the Phred quality use if the `-use-qual` flag was used)
* The entropy of the assignment

Note that the node identifiers in the `estmult.*` files are equal to (the BCALM node identifier + 1). In other words, we number nodes 1, 2, 3, ... whereas BCALM start numbering from 0. A negative node identifier refers to the reverse-complement node.


### 3. [OPTIONAL] Determining phased SNP calls with respect to a reference genome

The contigs should be long enough to be used with a typical variant calling tool to obtain strain-aware variant calls.
One posibility to obtain strain-aware variant calls (with respect to `reference.fasta`) is with the following procedure:

```
for ((i=1; i <= $nStrains; i++))
do
 idx=$((i-1))
 echo ">${idx}" >> $varfile
 PATH/TO/minimap2/minimap2 -cx sr -t4 --cs ../reference.fasta strain$i.phased.fasta | sort -k6,6 -k8,8n > strain${i}.phased.sort.paf
 PATH/TO/k8 PATH/TO/minimap2/misc/paftools.js call  -l 100 -L 100 -f ../reference.fasta -s "strain$i" "strain$i.phased.sort.paf" >> $varfile
done
```
## Advanced user information 

### The HaploDetox Pipeline in more detail
The HaploDetox pipeline consists of 4 stages. Each stage outputs a number of intermediate files. When re-running HaploDetox (e.g. using different parameter settings), the existence of these intermediate files is checked. If they exist, stages 1 or 2 may be skipped. To force re-running stages 1 or 2, simply remove the appropriate files: `rm *.st1` for stage 1,  `rm *.st2` for stage 2. We explain the stages in more detail below, and we list HaploDetox settings that are relevant for these stages.

#### Stage 1
The input file `reads.unitigs.fa` produced by BCALM 2 is read. Even though BCALM 2 provides average k-mer counts for each node (unitig), it does not provide (k+1)-mer counts for the arcs. Also, when using the `-use-qual` flag, HaploDetox needs to weigh the k-mer counts using the Phred quality scores. Note that the Phred score ASCII base value can be set using `-phred-base <value>` (default=33).
So in either case, after reading the de Bruijn graph into memory, HaploDetox needs to stream through all reads in order to establish both the appropriate node and arc counts. The results are written to `reads.unitigs.fa.st1`, a binary file containing the de Bruijn graph annotated with q- or k-mer coverage. This file is a prerequisite for stage 2.

If the file `strains.mf` containing a list of ground-truth fasta references is present in current working directory, HaploDetox will also compute the true node and edge multiplicities. These will be stored in files `truestrains.node` and `truestrains.edge`, respectively.

In principle, one should only re-run stage 1 when switching between k-mer or q-mer counts.

#### Stage 2
Model fit and graph cleaning

In stage 2, spurious nodes and arcs are removed from the de Bruijn graph. This is done iteratively. Based on an initial model fit to the k-mer spectrum, a maximum coverage value is selected below which a node or arc is considered as a low coverage node/arc `Cl`. Low coverage nodes and arcs are considered for possible removal. 
In a first pass, low coverage nodes/arcs are only removed when the conservation of flow property holds in all nodes in a subgraph around that node/arc (subgraph size set by `-crf-nb-size`). We split this pass up into several iterations such that the maximum considered coverage value is increased at each iteration until we reach `Cl`. After each iteration selected nodes/arcs are removed and the de Bruijn graph is contracted further where possible. In a second pass, the CRF of a subgraph is constructed for each (remaining) low coverage node/arc and the low coverage node/arc is only removed when its CRF-based multiplicity estimate is 0. This pass is again split up in sevaral iterations and selected nodes/arcs are removed and the de Bruijn graph contracted after each iteration.

 * `-no-correct` When the `-no-correct` flag is set, HaploDetox will skip graph correction. It is advised to correct the graph, however, as this improves the speed of the subsequent steps in HaploDetox significantly. It also improves accuracy when the number of variants is fairly low, but has the risk of removing very low coverage variants together with the sequencing errors. As a final advantage, the resulting contigs in the strain-aware output will be longer.

After graph correction, final mixture models are fitted that represent the k-mer coverage of nodes and arcs. We use Expectation-Maximisation with method of moment estimators and multiplicity estimates based on a Conditional Random Field.

If the estimated number of strains is not passed to HaploDetox, the process of graph-cleaning and model fitting is repeated under the assumptions of number of strains = 2, 3, and 4. Several model fit scores are output, based on which the most appropriate number of strains can be selected. We use the BIC score to select the number of strains, but other scores are output and you can rerun stage 3 and 4 with the strain number of your choice. Stage 3 requires `model.node.st2`, `model.edge.st2` and `reads.unitigs.fa.st2`. If you want to use another number of strains than was selected by the BIC score, simply copy the files `xstrains.model.node.st2`, `xstrains.model.edge.st2` and `xstrains.reads.unitigs.fa.st2` to the required file names for the number of strains x of your choice.

Settings that influence the mixture model fit are:

 * `-em-no-eqweight` We enforce approximately equal weights for components of the model that represent a multiplicity in which the same number of strains are present (e.g., for three strains (0,0,1), (0,1,0) and (1,0,0) will get similar weights, and (1,1,0), (1,0,1), (0,1,1) will also be enforced to have similar weights. This often leads to better convergence of model fitting, but can be disabled with the `-em-no-eqweight` option.

 * `-em-fix-zero` When this option is set, the error distribution will not be retrained after graph cleaning. Instead, it will be kept fixed to the error-component of the initial mixture model. Setting this option leads to worse performance of model selection scores in selecting the correct number of strains, but it might be useful when there is a strain with fairly low read support. It then avoids that the new error-component tries to include the k-mers belonging to this low read support strain.

 * `-mm-coverage`, `-mm-err-cov`, `-mm-odf`, `-mm-frac` With these options the initialisation of the k-mer spectrum EM can be set manually: use these to initialise respectively the avg. coverage $\lambda$ of unique $k$-mers common to all strains, the avg. error coverage, the negative binomial overdispersion factor and the strain-fractions $0<\pi_s<1$ such that the average coverage of strain $s$ is $\pi_s\lambda$. For setting `-mm-frac` make sure $\sum_s \pi_s = 1$ and pass the values as $\pi_1;\pi_2;...;\pi_S$.

Settings that influence (inference with) the CRF are:

 * `-no-approx-inf` By default, we use approximate inference (loopy belief propagation) on CRFs that represent the complete de Bruijn graph for the EM computations and in stage 3. We only use exact inference on subgraph-based CRFs during error correction. With this setting, approximate inference can be disabled such that exact inference on subgraph-based CRFs (set neighbourhood size with `-crf-nb-size`) is used for all multiplicity computations.

 * `-crf-nb-size` This option specifies the size of the neighbourhood subgraph around the node/arc whose multiplicity is inferred with exact inference (variable elimination). Larger neighbourhood sizes result in a longer runtime, but higher accuracy. The default value is 3. In most cases, using a small neighborhood size is enough. However, if regions of extremely high or low coverage span large regions of the de Bruijn graph, it might take a large neighbourhood size to obtain a CRF that can infer the correct multiplicities. If this results in too long runtimes, approximate inference can be used. 

 * `-crf-max-fact` This sets the maximum allowed intermediary factor size for variable elimination (exact inference). When you have enough RAM avaible and you want to run exact inference with a large neighbourhood size (this will increase runtime!) you might consider setting this to a higher value (default: $1e^{6}$).

 * `-crf-flow` This sets the value assigned to multiplicity combinations that adhere to conservation of flow of multiplicity within a strain. Other multiplicity combinations get a value of 1 (default: $1e^7$)

Settings that influence EM are:

 * `-em-max-iter` maximum number of EM-iterations (default: 25)

 * `-em-conv-eps` maximum relative change in log-likelihood to assume convergence (only used with approximate inference) (default: $1e^{-6}$)

 * `-em-max-change` maximum allowable proportion of training nodes/arcs that are changed between iterations before convergence is assumed (default: $1e^{-3}$).

 * `-em-train-size` number of nodes and arcs used for model training (set to -1 to train on all nodes and arcs) (default: $1e^4$).


#### Stage 3
Multiplicity determination of all nodes and arcs given the chosen number of strains using CRF.

#### Stage 4
Graph phasing and contig output.

#### Other option settings

 * `-num-threads` to specify the number of threads to use when using exact inference on subgraphs. By default, the number of threads equals the number of CPU cores or twice that number for CPUs that support hyperthreading.
 * `-help` to display the help page.
 * `-version` to display the version number.


### Checking the fit of the model produced in stage 2

If you want to check the fit to the histogram of the model that HaploDetox trains in stage 2, this can be done visually using [Gnuplot](http://www.gnuplot.info/).
Intermediate files `xstrains.haplotypeNBmix.gnuplot`, `xstrains.haplotypeNBmix.dat`, can be used to plot the histogram and the initial mixture model fit under the assumption of x strains. `xstrains.nodes.gnuplot` (resp. `xstrains.edges.gnuplot`) and `xstrains.nodes.dat` (resp `xstrains.edges.dat`) can be used to plot the histogram overlaid with the estimated mixture model under the assumption of x strains, after graph cleaning and re-training with EM and CRF-based multiplicity assignments. 

Just run
```
gnuplot prefix.gnuplot
```
This will produce a pdf output of the plot.
HaploDetox auto-determines a good scope for the plot, but you can manually change this in the `.gnuplot` files if need be.

#### What to do when the fit seems bad

We have determined what we think are the best default settings for the parameters HaploDetox uses, based on tests on several datasets. However, if your visual inspection of the model fit seems bad, you can try to manually alter some of these parameter values (see the list of settings in the description of stage 2)

##### Check if convergence was obtained

By default, the maximum number of EM-iterations is 25. You will be notified by HaploDetox when model fitting finished with this many iterations. If you see that there was still a relatively large change in parameter values between iteration 24 and 25, it can be beneficial to rerun stage 2 with a higher number of EM-iterations. Option `-em-max-iter [iters]` sets the maximum number of EM-iterations.

 ### Visualizing parts of the de Bruijn graph

You can visualize the initial de Bruijn graph using [Cytoscape](https://cytoscape.org/)

The graph is provided as two files `Cytograph.full.nodes` and `Cytograph.full.arcs`. Within Cytoscape you first import the `.arcs` file with `File > Import > Network from file...` and then you import the `.nodes` file with `File > Import > Table from file...`

You can also rerun HaploDetox with the option `-vis-subgraph` (to which you pass a nodeID). This replaces stage 3 and 4 with a routine that extracts a subgraph for a neighbourhood of size `-crf-nb-size`, computes all multiplicities within this subgraph, and exports files for Cytoscape visualisation (named `cytgraph[nodeID]nb[nb-size].nodes` and `cytgraph[nodeID]nb[nb-size].arcs`. This option requires the output files `*.st1` and `*.st2` from stage 1 and 2 to be present.
