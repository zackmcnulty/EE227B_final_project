# EE227B (Convex Optimization) : Final Project

With the growing demand for data capacity in communication, the emerging massive/large-scale multiple-input multiple-out (LS-MIMO) transmission techniques are employed in the forthcoming fifth generation (5G) communication systems for higher spectral-efficiency and energy efficiency. The challenges of accurate and efficient detection in the receiver side have become a popular topic in wireless communications and convex optimization has been used widely in design and analysis of communication systems and signal processing algorithms. Adopting semidefinite programming (SDP) relaxation to MIMO detection problem is one such usage. Compared to other existing detection algorithm which have exponential complexity such as sphere decoder,  SDP relaxation provides an effective polynomial-time worst-case computational complexity and offers an excellent performance-complexity tradeoff in practical SNR ranges. However, there are some constraints in terms of the modulation scheme and the channel model. For SDP relaxation to work properly, the channel would needed to be real-valued while it is natural to have a complex channel. On the other hand, for channel with memory, one-shot detection is not applicable and some extension of the matrices could be done to help solve the sequence detection with hardware overhead and computation penalty.

In this project, we first briefly introduce the MIMO system and define the MIMO detection problem. From there, we review several popular pre-existing algorithms used for MIMO detector, including maximum likelihood (ML), zero-forcing (ZF), minimum mean-square error(MMSE), and an SDP relaxation of the ML problem. Based on the SDP reformulation, we further propose a new heuristic algorithm for handling the nonconvex rank one constraint. Furthermore, we explore some possible robustifications of the maximum-likelihood MIMO detection problem that results from the inability to exactly measure some of the MIMO system parameters. 


### Repository Contents

###### robust_ml

Code used for generating figures for the robustification of the maximum likelihood (ML) problem under L1 and L2 norm ball uncertainty (sections 4.1.1 and 7.3 of the report).

###### existing_algorithms

Code used for reproducing results from pre-existing algorithms for the MIMO detection problem (sections 3 and 7.1 of the report).
