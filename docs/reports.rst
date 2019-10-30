Reports
~~~~~~~

Summary Report
^^^^^^^^^^^^^^
Large .tsv file containing detailed information for each cluster detected by centreseq.
Displays representative sequence label, number of members, and which sequence labels belong to the cluster.

Core Gene Count Report
^^^^^^^^^^^^^^^^^^^^^^
Contains general metrics on the # of core genes detected among 100% of samples, >=90% of samples, and >50% of samples.

Roary Gene Count Report
^^^^^^^^^^^^^^^^^^^^^^^
Gene count report represented in the style of Roary's output.

Pairwise Report
^^^^^^^^^^^^^^^
Large .tsv file which stores pairwise information between samples on the number of matching genes and non-matching
genes among the intersection of all genes existing between the two. Used to generate the network visualization.

Network Chart
^^^^^^^^^^^^^
Interactive visualization generated from the pairwise report file.

.. image:: images/network.png
  :width: 300
  :alt: Example network chart

Rarefaction Curve
^^^^^^^^^^^^^^^^^
Visualization generated showing the # of 'core' genes vs. the # of 'pan' genes with increasing numbers of
sampled genomes.

.. image:: images/rarefaction_curve.png
  :width: 600
  :alt: Example rarefaction curve

Looping over a range from 1..n, the following process is executed,
    1) n samples are randomly selected
    2) Number of core genes shared between the subset calculated
    3) Total number of genes existing among the subset is calculated

By default, this entire process is repeated 5 times to reduce variance from the random sampling.
