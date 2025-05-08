# SfQSM
SfQSM is a three-dimensional tree modeling method based on skeleton graph optimization and fractal self-similarity. This method consists of three main steps:ⅰ) Skeleton points self-adjusting based on geometric features, ⅱ) Edge weight definition for the tree graph, and ⅲ) Fractal self-similarity optimization for individual tree modeling. 

More detailed information about SfQSM can be found in the article "Self-adaptive Individual Tree Modeling Based on Skeleton Graph Optimization and Fractal Self-similarity." 

If you use this code, please remember to cite this paper.

# Code structure

SfQSM.m : The main function.

Self-adjustment : This directory contains the fundamental function codes for adjusting skeleton points based on geometric features. 

tree_graph : This directory includes the essential function codes for constructing individual tree graphs with a new edge weight definition. 

TreeModeling : This directory encompasses the fundamental function codes for modeling individual trees based on Fractal self-similarity optimization. 

toolbox : This directory includes the basic function.

skeleton : This directory contains the skeleton extraction code based on Laplace contraction, originally developed by Teacher Cao Junjie from Dalian University of Technology. Our method has been modified and built upon this technique.

Data : This directory holds the test data in '.txt' format. 

Result : This directory stores the modeling results of the test data. 

plotting : This directory includes the basic function codes for visualization from TreeQSM.

The SfQSM is programmed in Matlab R2022b.
