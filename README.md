# SPaRTAN
**SPaRTAN (Single-cell Proteomic and RNA based Transcription factor Activity Network)** is a computational method to link surface proteins to transcription factors (TFs)  in which we exploit cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq) datasets with cis-regulatory information.

<img src="https://github.com/osmanbeyoglulab/PyAffreg/blob/master/data/diagram.png" width="500">

Our model starts from a readout of surface protein activities as protein expression, surface activity then passes signals along the pathway and converges on a set of TFs, whose activity, in turn, regulates the mRNA expression levels of TF target-genes.
We observed gene expression(Y) as the interaction between two profiles, protein expression(P) and TF against target genes matrix(D). We utilized those relationships and affinity regression to establish an interaction matrix between surface proteins and TFs (W) and then further predicted  target gene expression.

SPaRTAN can be simplified as solving the equation

**DWP<sup>T</sup>=Y**

in which D, P, Y are the data input and W is the matrix we need to solve.

We first conduct dimension reduction by multiplying Y<sup>T</sup> on both sides of the equation, then use SVD (singular vector decomposition) to further reduce the dimension of matrices. After applying a series transformation, we convert this bilinear problem into linear regression, in which 4 parameters: pectrumA, spectrumB, rsL2 and lambda need to be tuned based on the user inputs, where spectrumA and spectrumB are related with SVD; rsL2 and lambda are related with linear regression.

There are currently implementations of SPaRTAN in Matlab, and in Python. For more details and installation instructions on running SPaRTAN, please see the tutorials
* [Run SPaRTan in Python](https://github.com/osmanbeyoglulab/SPaTRAN2/tree/main/SPaRTAN_python)
* Run SPaRTan in Matlab
