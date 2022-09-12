# tDigestForEcon
Using t-Digest to calculate quantiles of agent distribution (codes for paper)

These are the matlab codes to accompany paper: Robert Kirkby - Computing Quantiles for the Agent Distribution using t-Digests

Two examples:
 - tDigest.m: calculates the quantiles both exactly and using t-Digests for a few different randomly generated samples/distributions.
 - LifeCycleModel_tDigest.m: solves a life-cycle model and calculates the quatiles of functions of agent distribtuion both exactly and using t-Digest

The implementation of the t-Digests is done by the two functions:
 - createDigest(): creates t-Digest from a distribution (of values and weights)
 - mergeDigest(): create t-Digest by merging together a few t-Digests
Both of these two functions are just copies of those included in VFI Toolkit

To solve the life-cycle model you will also need a copy of VFI Toolkit: vfitoolkit.com
(The version of VFI Toolkit at end of Sept 2022 is the one that was used for the paper.) 
