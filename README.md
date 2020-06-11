# Gene-Annotation-and-Prediction
This project is about predicting the location of a gene in the human genome based off experimental data. 

![DNA Strand](https://github.com/jaotheboss/Gene-Annotation-and-Prediction/blob/master/dna.jpg)

## Objective:
1. Predict the location of a gene in the human genome from a BAM file.

2. Compare the predictions with the TxDb file to evaluate the predictions and the model.

## Context:
Gene annotation is incredibly difficult as there are multiple ways to detect the start and end of the gene within a chromosome. However, for someone with no background in life sciences whatsoever, but is still keen on practicing data in this domain in some way, this simple project would suffice.

Within the cell of every human cell, there are 23 pairs of chromosomes (there exists 24 kinds of chromosomes), with 22 pairs being of the same kind and the 23rd pair (the sex chromosome) being either a pair of X chromosomes or a X-Y pair. These chromosomes are made out of DNA (an acid molecule). Sections of these DNA are called genes, and the collection of all unique genes is called a genome; hence, the human genome is just a collection of genes that make up a human being. Genes are essentially instructions for our cells to form/create whatever we are now. Ranging from the size of your nose, color of your eyes, to the color of your skin, these traits are instructed by genes within your cell. 

There is a difficulty in being able to detect which parts of the DNA is considered a gene, due to the sheer complexity of the subject. Hence, even readings from a specialized piece of equipment could vary largely. Moreover, that's not to say that the piece of equipment used is 100% perfect as well; the equipment itself may have accuracy errors.

Therefore, for this project, i will be skipping the reading of the DNA sequences to predict where the gene starts and ends, but instead i will be processing the reading results from the equipment and attempt to predict the start and end locations from there. The reading results consists of the start and end indexes of where the equipment decides there to be the start and end of a gene. Hence, the equipment can have multiple readings for 1 gene. 

## Dataset:
In this project, i will essentially be working with 2 types of data: BAM and TxDb.(Both of which are included in the zip file of this repo)

1. BAM

DNA sequences are usually saved in a file format called Fastq. However, this file format is only useful for ascertaining the quality of the DNA sequence reading itself. Hence, in order to work on the data, the Fastq file has to be transformed into something more standard. This Fastq file is thus compared with a reference genome (holds information about the genes we're interestd in; like an answer sheet), labelled at their start and end, and then aligned to eventually form a SAM format file. The BAM and SAM file format is essentially the same, just that BAM is binary and more data analytically friendly. 

Hence, the BAM file looks like this:

seqnames | strand | start | end | width
-|-|-|-|-
chr1|+|1|100|100

Note: Further explanation of the column names. **seqnames** = indicates which chromosome this gene is from, **strand** = since DNA is a double helix, the + and - indicate which side is being used as a gene, **start, end** = start and end of the gene, **width** = length of the gene

2. TxDb

The TxDb file is essentially a container of information regarding a certain genome. This can be thought of as the answer booklet of the genes for that genome. This container works the same way as an SQLQuery, where we can query different information from it. We can extract where the genes are located, which parts are the exons, introns and etc.

## Method

### Model:
The model in the Python script first splits the genome into their respective chromosomes and strands. The model then works on this smaller dataset one at a time before returning a file that’s a compilation of the results from all these smaller datasets. 

For each of these datasets, the start and end of the reads are sorted. After which, all the starting indexes were clustered up iteratively, only forming a cluster if the next starting index was statistically improbable to be within said cluster. This was done through using a normal distribution of a moving average that’s made up with those currently in said cluster, and a fixed variance. 

The number of clusters formed are the number of genes that the model predicts to be present. The fixed variance and the threshold for improbability (p_value) can be chosen before the execution of the script. 

### Evaluation: 
This model and its result is evaluated by comparing every gene prediction with the gene indexes from the TxDb file, following a certain threshold. This threshold can also be chosen before the execution of the evaluation. The threshold that’s currently set is 20692 as that’s the average length of a gene (based on the data in the TxDb file). This implies that as long as the gene prediction is plus minus 20692 index away from the actual, they would be considered accurate. Hence, if we were to look for a more accurate result, we would have to reduce this threshold. 
