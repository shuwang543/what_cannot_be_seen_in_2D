

# Supplementary Information to “What Cannot Be Seen Correctly in 2D Visualizations of Single-Cell ‘Omics Data?” by Shu Wang, Eduardo D Sontag, Douglas A Lauffenburger.

## Section S1: Curvature affects geodesics and their representations

Generically, it is impossible to map an intrinsically curved manifold to an intrinsically flat space such as a 2D plane or other Euclidean space in such a way that the manifold’s geodesics can all be represented as straight lines in the map. This informal statement is the essence of Beltrami’s Theorem (1) and the Killing-Hopf Theorem (2), which are two of many theorems characterizing how curved spaces lack equivalent representations in flat spaces. The exceptions to the generic case only occur when a manifold has a constant value of curvature at every point, in addition to other requirements. For example, the spherical surface of the Earth is a standard example of a manifold with constant curvature at every point, although it fails the other requirements (i.e. the sphere cannot be smoothly mapped to a plane due to discontinuities). However, it is technically possible to map a small patch of the Earth such that local geodesics can still be qualitatively represented as straight lines (e.g. the Gnomonic projection), although at the expense of distorting quantitative properties such as distances, areas and angles (e.g. the equator would be infinitely far away). Regardless, since we expect nonlinear data manifolds to be generic, intrinsically curved manifolds, any 2D map of it will fail to qualitatively represent the shortest paths between data points as straight lines.

Geodesics have been interpreted as paths of cell-type differentiation (3) or organismal evolution (4), and recent work has also used concepts of manifold Laplacians and diffusion distances to investigate differential gene expression (5) and to infer “pseudo-time” progression (6) in single-cell data. Other concepts from differential geometry have corresponding data-oriented computational tools that may also be useful for framing biological hypotheses: metrics (7), curvatures (8), the profile of diffusive processes (9), and other features can all behave in unintuitive fashions when the data manifold is nonlinear. 

Notably, when one considers how measurement noise affects data points sampled from a nonlinear manifold, it quickly becomes intractable to recover the original manifold compared to the simpler case of linear manifolds (10–12). For example, in the case of finding the principal curve (one of the nonlinear generalizations of the 1st Principal Component) from sampled data points, there is a statistical bias of the fitted principal curve relative to ground-truth, dependent on the ratio of noise to curvature (13). On a more optimistic note, if there is prior knowledge to believe that the original data manifold falls into certain special classes, there are corresponding transformations (e.g. as simple as the logarithm) and algorithms that can help in recovering a specific manifold model of a data set (24–26).

## Section S2: Discrete and combinatorial geometry applied to single-cell biology

Given four equidistant points, it is impossible to place them simultaneously onto the Euclidean plane while maintaining equidistance. Furthermore, the deviation from equidistance would not be infinitesimal, as shown in Fig 1B, and the deviation grows worse if one has more than four points that are meant to be equidistant (14). Also, if distance in the N-dimensional data space is not measured by Euclidean distance, but instead by other metrics, there are potentially more than N+1 equidistant points; even if data points are constrained exactly to a two-dimensional surface, if the surface is sufficiently curvy, there could still be k equidistant points, for arbitrarily large k, if one considers geodesic distances on this surface (15), and so despite the data distribution being “two-dimensional”, one could not construct a representation on the 2D Euclidean plane that preserves equidistance of those k>3 points.

Problems concerning the possibility of arranging points to fulfill various qualitative properties like distance inequalities are in the tradition of combinatorial or discrete geometry. Equidistant point arrangements can be considered a subset of problems like the Erdös distance problem (16), in which instead of determining if k points can be arranged to have a unique pairwise distance d (i.e. equidistant) one seeks point arrangements so that the set of distinct pairwise distance values {di} is as small as possible. Alternatively, the problem of how the maximum k depends on the metric of a space (e.g. the Chebychev, Minkowski, or Riemannian metrics that are common choices in single-call data analysis) is studied as a problem of equidistant sets (17, 18). Meanwhile, the issue of enumerating possible distance permutations is directly related to the classic combinatorial problem of cutting a cake into pieces with k slices (19): each piece is a polygonal region, analogous to the regions shown in Fig 2A cut out by dotted lines. Generalizations of these mathematical results can have value for applications as well: the capacity of different machine learning algorithms to classify objects into discrete categories is often quantified by their Vapnik-Chervonekis dimension (20, 21), which is calculated based on these same sets of combinatorial tools. Discrete geometry objects like simplices and polytopes have also been used to model single-cell phenotypes in different tissues (22, 23), the spatial composition of tissues in terms of micro-environmental niches (24), and evolutionary fitness landscapes with epistasis (25).

## Section S3: Topological data analysis of biological data 

Topological data analysis (TDA) aims to infer data properties that are robust to noise (26) by focusing on topological descriptors, which by definition can only change when discontinuous events such as cutting or gluing occur. Given a dataset, TDA offers many rigorous ways to define topology from data points, each with various pros and cons (27), which have recently been used to identify novel candidate cancer-associated genes in various tumor types (28), and for understanding the evolutionary relationships of species beyond tree topologies (29, 30). Similarly, the past two decades of hematopoiesis research has begun to challenge the canonical tree-model’s branching topology (31–33) of immune cell development, and TDA would provide natural tools for characterizing the inherent topology of hematopoiesis.

## References:
1. H. Busemann, B. B. Phadke, A general version of Beltrami’s theorem in the large. Pac. J. Math. 115, 299–315 (1984).
2. J. M. Lee, Introduction to Riemannian Manifolds (Springer, 2019).
3. G. La Manno, et al., RNA velocity of single cells. Nature 560, 494–498 (2018).
4. L. J. Billera, S. P. Holmes, K. Vogtmann, Geometry of the Space of Phylogenetic Trees. Adv. Appl. Math. 27, 733–767 (2001).
5. K. W. Govek, V. S. Yamajala, P. G. Camara, Clustering-independent analysis of genomic data using spectral simplicial theory. PLoS Comput. Biol. 15, e1007509 (2019).
6. X. Sun, J. Zhang, Q. Nie, Inferring latent temporal progression and regulatory networks from cross-sectional transcriptomic data of cancer samples. PLoS Comput. Biol. 17, e1008379 (2021).
7. S. ren Hauberg, O. Freifeld, M. Black, A Geometric take on Metric Learning in Advances in Neural Information Processing Systems, (Curran Associates, Inc., 2012).
8. D. Sritharan, S. Wang, S. Hormoz, Computing the Riemannian curvature of image patch and single-cell RNA sequencing data manifolds using extrinsic differential geometry. Proc. Natl. Acad. Sci. 118 (2021).
9. M. Belkin, P. Niyogi, Towards a theoretical foundation for Laplacian-based manifold methods. J. Comput. Syst. Sci. 74, 1289–1308 (2008).
10. P. Niyogi, S. Smale, S. Weinberger, Finding the Homology of Submanifolds with High Confidence from Random Samples. Discrete Comput. Geom. 39, 419–441 (2008).
11. C. Genovese, M. Perone-Pacifico, I. Verdinelli, L. Wasserman, Minimax Manifold Estimation. J. Mach. Learn. Res. 13, 1263–1291 (2012).
12. C. R. Genovese, M. Perone-Pacifico, I. Verdinelli, L. Wasserman, Nonparametric ridge estimation. Ann. Stat. 42, 1511–1545 (2014).
13. T. Hastie, W. Stuetzle, Principal Curves. J. Am. Stat. Assoc. 84, 502–516 (1989).
14. N. Linial, E. London, Y. Rabinovich, The geometry of graphs and some of its algorithmic applications. Combinatorica 15, 215–245 (1995).
15. R. K. Guy, An Olla-Podrida of Open Problems, Often Oddly Posed. Am. Math. Mon. 90, 196–200 (1983).
16. J. Garibaldi, A. Iosevich, S. Senger, The Erdös Distance Problem (American Mathematical Society, 2011).
17. C. M. Petty, Equilateral Sets in Minkowski Spaces. Proc. Am. Math. Soc. 29, 369–374 (1971).
18. J. Wilker, Equidistant sets and their connectivity properties in (1975) https:/doi.org/10.1090/S0002-9939-1975-0355791-0.
19. R. C. Buck, Partition of Space. Am. Math. Mon. 50, 541–544 (1943).
20. T. M. Cover, Geometrical and Statistical Properties of Systems of Linear Inequalities with Applications in Pattern Recognition. IEEE Trans. Electron. Comput. EC-14, 326–334 (1965).
21. E. D. Sontag, “VC Dimension of Neural Networks” in Neural Networks and Machine Learning, (Springer, 1998), pp. 69–95.
22. Y. Korem, et al., Geometry of the Gene Expression Space of Individual Cells. PLOS Comput. Biol. 11, e1004224 (2015).
23. M. Adler, Y. Korem Kohanim, A. Tendler, A. Mayo, U. Alon, Continuum of Gene-Expression Profiles Provides Spatial Division of Labor within a Differentiated Cell Type. Cell Syst. 8, 43-52.e5 (2019).
24. A. E. Marrahi, F. Lipreri, D. Alber, J. Hausser, Four tumor micro-environmental niches explain a continuum of inter-patient variation in the macroscopic cellular composition of breast tumors. 2022.03.04.482793 (2022).
25. H. Eble, M. Joswig, L. Lamberti, W. B. Ludington, Cluster partitions and fitness landscapes of the Drosophila fly microbiome. J. Math. Biol. 79, 861–899 (2019).
26. G. Carlsson, Topology and data. Bull. Am. Math. Soc. 46, 255–308 (2009).
27. N. Otter, M. A. Porter, U. Tillmann, P. Grindrod, H. A. Harrington, A roadmap for the computation of persistent homology. EPJ Data Sci. 6, 1–38 (2017).
28. R. Rabadán, et al., Identification of relevant genetic alterations in cancer using topological data analysis. Nat. Commun. 11, 3808 (2020).
29. J. M. Chan, G. Carlsson, R. Rabadan, Topology of viral evolution. Proc. Natl. Acad. Sci. 110, 18566–18571 (2013).
30. P. G. Cámara, A. J. Levine, R. Rabadán, Inference of Ancestral Recombination Graphs through Topological Data Analysis. PLOS Comput. Biol. 12, e1005071 (2016).
31. R. Ceredig, A. G. Rolink, G. Brown, Models of haematopoiesis: seeing the wood for the trees. Nat. Rev. Immunol. 9, 293–300 (2009).
32. S. Watcham, I. Kucinski, B. Gottgens, New insights into hematopoietic differentiation landscapes from single-cell RNA sequencing. Blood 133, 1415–1426 (2019).
33. S. Haas, A. Trumpp, M. D. Milsom, Causes and Consequences of Hematopoietic Stem Cell Heterogeneity. Cell Stem Cell 22, 627–638 (2018).

