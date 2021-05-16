# SPaRTAN
Signal-regulated transcription factors (TFs) underlie key developmental or differentiation transitions and activation states of cells (e.g., within the immune system). To date, single-cell genomics datasets have not been used to link cell-surface receptors to TFs of individual cells. To address this need, we developed a machine learning approach called **SPaRTAN (Single-cell Proteomic and RNA based Transcription factor Activity Network)** that integrates parallel surface protein measurements with mRNA expression data in cells based on CITE-seq (cellular indexing of transcriptomes and epitopes by sequencing) data with regulatory genomics resources. 

<img src="https://github.com/osmanbeyoglulab/SPaRTAN/data/diagram.png" width="500">

Briefly, our model views expression of surface proteins as a proxy of their activities; signaling emanating from these proteins converges on particular TFs, whose activities, in turn, regulate the expression of their target genes. Specifically, we use a regularized bilinear regression algorithm called affinity regression (AR), a general statistical framework for any problem where the observed data can be explained as interactions (**W**) between two kinds of inputs, to establish an interaction matrix between surface proteins (**P**) and TFs (**D**)  that predicts target gene expression (**Y**). 

Since the model captures statistical relationships between surface proteins, TFs, and gene expression. We can use the trained interaction matrix to obtain different views of a CITE-seq data set; e.g., to predict TF activity from a cell's surface protein expression profile or to predict surface protein expression from a cell’s gene expression profile.  Intuitively, information flows down from observed surface protein levels through the learned interaction matrix to infer TF activities and observed mRNA expression levels or propagates up through the TF target-gene edges and interaction network to infer surface protein expression. 

There are currently implementations of SPaRTAN in Matlab, and in Python. For more details and installation instructions on running SPaRTAN, please see the tutorials
* [Run SPaRTAN in Python](https://github.com/osmanbeyoglulab/SPaTRAN2/tree/main/SPaRTAN_python)
* Run SPaRTAN in Matlab
