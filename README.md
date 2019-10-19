# Finite-graph-fractal-sets


## The big idea of what this code can do
This repository contains python code for creating fractal sets from labellings on a finite directed graph.
We use the tribonacci graph with the following labelling on the edges.

<img src="graph.png" alt="alt text" width="250" height="250">

A user can decide on a labelling (a,b,c) where a,b,c are real numbers, to each of the edges of the graph and the code will output a fractal-like set in the xy-plane.  

Here are a few examples of the sets that the code produces.  

![edge4](edgelabeling1.png)


![edge4](edgelabeling2.png)


We can combine the above two labellings to create the Rauzy fractal!
![edge4](edgelabelingrauzy.png)

Just as another example:

![edge4](edgelabelingextra.png)

##The code needs numpy, matplotlib,itertools

## The functions involved
This code contains five main functions.

The first functions calculate paths edges of a desired length on the finite directed graph, forward, backward and begining at all three vertices.  The second type of function relabeles all of the paths according to the chosen labelling.  The third function applies the specfic map from my research to each labelled path.  The fourth functions project the path either on the contracting plane of the hyperbolic matrix or the expanding line.  The last main function projects points on the contracting plane to an xy plane for easier viewing.  

##Mathematical questions

The mathematics behind this code involves Markov partitions of hyperbolic toral automorphisms.  There are a lot of mathematical questions to explore with this code.  Which labelings produce sets which are the closure of their interiors? Which labelings produce totally disconnected sets?  Which labelings give non-overlapping sets? 
