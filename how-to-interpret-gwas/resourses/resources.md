---
description: '来自知乎橙子牛奶糖 https://zhuanlan.zhihu.com/p/90414014'
---

# 全基因组关联分析学习资料（GWAS tutorial）

## 前言

很多人问我有没有关于全基因组关联分析（GWAS）原理的书籍或者文章推荐。

其实我个人觉得，做这个分析，先从跑流程开始，再去看原理。

为什么这么说呢，因为对于初学者来说，跑流程就像一个大黑洞，学习原理就像一个小黑洞。

很多人花了好几个月的时间在看原理，一旦丢给他数据去分析，依旧束手无策。

不会跑流程，内心依旧会很恐慌。就像从来没有入门一样。

所以，我的建议是咱们先不去管原理，直接从分析入手。

等把数据跑出来了，整个流程的技能点满了，再去看看它的原理。

## 入门：学习GWAS的在线网站：

对于没有编程基础的人来说，建议先从一个在线的网站走一遍GWAS流程。

这样就能知道完成GWAS需要多少个步骤，心里大概有个底。

[easygwas](https://link.zhihu.com/?target=https%3A//easygwas.ethz.ch/)网站提供了公共数据，可以直接开始分析GWAS。整个流程按照网站提示，很简单。

网址：[https://easygwas.ethz.ch/](https://link.zhihu.com/?target=https%3A//easygwas.ethz.ch/)

## 进阶备选1：在linux下学习GWAS的实操数据

由于我们最终还是需要拿着自己的数据完成GWAS分析，不必避免的需要一定的编程基础。

在线网站只是一个提供理解GWAS流程的网站，因此，我们还是需要在linux系统下拿一些数据练练手。学会最基本的命令行。

在这里，我推荐一个提供linux下学习GWAS的教程：[GWA\_tutorial](https://link.zhihu.com/?target=https%3A//github.com/MareesAT/GWA_tutorial/).

网址：[https://github.com/MareesAT/GWA\_tutorial/](https://link.zhihu.com/?target=https%3A//github.com/MareesAT/GWA_tutorial/)

网站分为四个教程：1）GWAS的数据QC; 2\)处理群体分层； 3）关联分析（GWAS）; 4\)多基因风险得分分析（Polygenic risk score analyses）

**示例数据都有了，就等你自己上手了。**

**我敢保证，当你能完整的跑完这个流程的时候，你对GWAS的理解少说也有70% ，下一个在群里帮我解答问题的大神就是你了（申请进群方式见公众号菜单栏）。**

## 进阶备选2：使用R语言做GWAS分析

有些人对R语言可能比较熟悉，这里提供了一个用R语言分析GWAS的流程。

该流程有：GWAS的QC，PCA分析，Manhattan图，QQ图，候选位点的功能分析

感兴趣的看这个：[Genome-wide association studies in R](https://link.zhihu.com/?target=https%3A//www.r-bloggers.com/genome-wide-association-studies-in-r/)

网址：[https://www.r-bloggers.com/genome-wide-association-studies-in-r/](https://link.zhihu.com/?target=https%3A//www.r-bloggers.com/genome-wide-association-studies-in-r/)

## 进阶备选3

0 原理

[啊，全基因组关联分析（GWAS）的计算原理，了解一下？](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483991%26idx%3D1%26sn%3Dc69e4db6124d6cafec175d529a05aa46%26chksm%3Dce2d6c37f95ae521a84bc374e30ce2c6a8fa57d64156f79af542d5de123ea0840c4049b82b88%26scene%3D21%23wechat_redirect)

1 分析流程

[GWAS分析基本流程及分析思路](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483663%26idx%3D1%26sn%3Dacfacdf0a0ee6df2c003875c0db06476%26chksm%3Dce2d6f6ff95ae67902343fd81768fa0e3cf48949d8ea01e1c3d50f4e09763e779906fda8de54%26scene%3D21%23wechat_redirect)

2 数据处理

2.1 数据质量过滤

[GWAS基因芯片数据预处理：质量控制（quality control）](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483780%26idx%3D1%26sn%3Dba5341ab9caa99a9a62b6810bcb48262%26chksm%3Dce2d6fe4f95ae6f2bbf7899a66acd9a0ad075d4db3159a7f9ab53e8668503c813b79c74a842c%26scene%3D21%23wechat_redirect)

2.2 正负链翻转（stand flip）

[数据合并，踩不完的坑](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483963%26idx%3D1%26sn%3Df1506df98c6578921ec68f522694295a%26chksm%3Dce2d6c5bf95ae54d6585498eca1c661eecbf58a30e757257f83b1753373471a19ce0e2eb501a%26scene%3D21%23wechat_redirect)

2.3 基因型数据填补（imputation）

[soga，网页版的基因型填充可以这么做（genotype imputation）](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483875%26idx%3D1%26sn%3D26e5656e521d6cb392291dd38d89de63%26chksm%3Dce2d6f83f95ae695458e80f5f0407ee596d69f0eae5e2242e9f848170e6136bb0a37fc4b0739%26scene%3D21%23wechat_redirect)

2.4 群体分层校正

[GWAS群体分层 \(Population stratification\)：利用plink对基因型进行PCA](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483767%26idx%3D1%26sn%3Dbc45f24c8d0e2369462fd9531be34eb6%26chksm%3Dce2d6f17f95ae601f4ca8a4d108654d7bb6324519a9020c4bf879ddce05ba82b63042adff725%26scene%3D21%23wechat_redirect)

[EIGENSTRAT除了用来计算PCA，还可以干嘛](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483944%26idx%3D1%26sn%3Df38f2b832dbbbcb479b84de5336c8d76%26chksm%3Dce2d6c48f95ae55e09cd720682b0b563686719380765bbf809be61fce73fd5f69ad4b4d51567%26scene%3D21%23wechat_redirect)

[群体遗传分析分层校正，该选用多少个PCA?](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483869%26idx%3D1%26sn%3D7e6797c8c7575a86abe239e3c26cfd63%26chksm%3Dce2d6fbdf95ae6ab46c8893737c2d5e850810153681be59315b6faa9ae238c25b932438403ac%26scene%3D21%23wechat_redirect)

3 关联分析

[GWAS: 曼哈顿图，QQ plot 图，膨胀系数（ manhattan、Genomic Inflation Factor）](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483759%26idx%3D1%26sn%3D4069ccf557f2efcdba3e5b19b0d20393%26chksm%3Dce2d6f0ff95ae619cb695d6f96b77ede84306f1838447385a71dd1cb27ba55c8e08eb9591ce9%26scene%3D21%23wechat_redirect)

4 meta分析

[只用一行命令，就可以学会全基因组关联分析\(GWAS\)的meta分析](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483888%26idx%3D1%26sn%3D79b8a2b07505c14fa4ded1d8ca2398a7%26chksm%3Dce2d6f90f95ae6861d9f8a74ee71c8291a0e6fe5e49fc4e125079224260aadfe991d8ebae4d8%26scene%3D21%23wechat_redirect)

5 条件分析

[GWAS条件分析（conditional analysis）：作用，步骤，结果解读](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483742%26idx%3D1%26sn%3D60fa1f3abea2e29898dd0be603195240%26chksm%3Dce2d6f3ef95ae628cd50d6dd52bbcf5c43ca0d7077b518b8cf8cd7fc938adaa12b77e7f17941%26scene%3D21%23wechat_redirect)

6 基因多效性

[啊啊救救我，为何我的QQ图那么飘（全基因组关联分析）](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483913%26idx%3D1%26sn%3Dbcde50e5d75da376a37a4fe885fd43cb%26chksm%3Dce2d6c69f95ae57f08c960a897bc1420e2ab1d66ec32abfc407605feb6b1f46bbda778ba39c2%26scene%3D21%23wechat_redirect)

7 GWAS后续分析

[SNP在世界地图上的频率分布](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483712%26idx%3D1%26sn%3Dc408554d411744bd64e850c648975d27%26chksm%3Dce2d6f20f95ae636de39e0b00f51bb32c0a22c44cd79e811267c7eb845c46e4b65d51b1a2804%26scene%3D21%23wechat_redirect)

[GWAS：拒绝假阳性之case和control数量比例严重失衡的解决方案（SAIGE模型的应用）](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483823%26idx%3D1%26sn%3D3aab057f961890546cb7d263f19d6da4%26chksm%3Dce2d6fcff95ae6d96e1526d7ba81f6d7671edd62b60bd48545114983483d5d108e384b0a0a5d%26scene%3D21%23wechat_redirect)

[GWAS后续分析：LocusZoom图的绘制](https://link.zhihu.com/?target=https%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483791%26idx%3D1%26sn%3D913f76042861b5b45ee9dc908c589d9b%26chksm%3Dce2d6feff95ae6f96007cca91136e69c1c0635654fbdef07f0a436779340d72dec0a5694ef5d%26scene%3D21%23wechat_redirect)

[利用GCTA工具计算复杂性状/特征（Complex Trait）的遗传相关性（genetic correlation）](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483682%26idx%3D1%26sn%3D0f0190a7f6a6d38db3e6eb5437e71ed4%26chksm%3Dce2d6f42f95ae654824bcbe7ccef0125ab8a083508740ee15d547d294b89dec5a98d40692625%26scene%3D21%23wechat_redirect)

[GWAS系列分析：多基因风险评分\(Polygenic Risk Score\)的计算](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483853%26idx%3D1%26sn%3D9f094b44447467d800fbbca01ca98b33%26chksm%3Dce2d6fadf95ae6bb997c2e2807914fc2ede4aaebf5f5037ddb51ed36d686d8c4b4ac5b5a6828%26scene%3D21%23wechat_redirect)

[有相关性就有因果关系吗，教你玩转孟德尔随机化分析](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483961%26idx%3D1%26sn%3Dee2c8afbc86361165fde3ffc7479f205%26chksm%3Dce2d6c59f95ae54fce8c46e438b27ef49c85b0fe88704144de4a4c5a2300da793f3cf2c9dc03%26scene%3D21%23wechat_redirect)

[LD SCore计算基因多效性、遗传度、遗传相关性](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483972%26idx%3D1%26sn%3D6abb84404d5c1e98ef496d9ee71c688e%26chksm%3Dce2d6c24f95ae532c5d12fec895763dafbf3c9afe6a8056b3c044621d5f7a9befd01f3941438%26scene%3D21%23wechat_redirect)

[查询、下载GWAS目录数据的R包\(gwasrapidd\)](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483978%26idx%3D1%26sn%3Db376a321362ae96d00b391656be61481%26chksm%3Dce2d6c2af95ae53c8cf03a5c266a3d682faff1db635939d502dbb667d1bf497361a4cb578e63%26scene%3D21%23wechat_redirect)

[DEPICT实现基因优化、gene set富集分析、组织富集分析（tissue enrichment）](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483985%26idx%3D1%26sn%3D782310edcec09c13ab5405f7ef06810a%26chksm%3Dce2d6c31f95ae52741330c1974fc7560a9e15861574c9cdd7b2152e68ea71471aabb205ece0b%26scene%3D21%23wechat_redirect)

[SNP功能注释网站合集](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247484028%26idx%3D1%26sn%3D3073a0a1d442889ffa5428881442335b%26chksm%3Dce2d6c1cf95ae50a1014ed6b0a0dd6177b0992e977502d4b8ba254cc417da25010be22a60452%26scene%3D21%23wechat_redirect)

8 相关文献阅读

[常见变异影响常见疾病，罕见变异影响罕见疾病？](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483880%26idx%3D1%26sn%3D5dc82bde1859d950f1ce288c9c51d021%26chksm%3Dce2d6f88f95ae69e0027152a3f5435bf9d7558e8dfa7040b7062453f4a921bb88e1fecbae2de%26scene%3D21%23wechat_redirect)

[什么！GWAS研究中case和control的比例是有讲究的？](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483809%26idx%3D1%26sn%3Dcbfda3d91f0453bb4720ff9f05431c02%26chksm%3Dce2d6fc1f95ae6d792a7d305a04bb9b9dc3d99ad110b4f4dfe6ed38dce2932147ea74a4effc4%26scene%3D21%23wechat_redirect)

[阿尔兹海默症和代谢指标在大规模全基因组数据的遗传共享研究](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483774%26idx%3D1%26sn%3Dbbfa8638f877c401ec26b6fec7e3fb6c%26chksm%3Dce2d6f1ef95ae608ce7e3aa007676bd4a012531f827563e8860c6f71d6a13e0239a91096090b%26scene%3D21%23wechat_redirect)

[GWAS文献解读：The stability of educational achievement](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483731%26idx%3D1%26sn%3D1f8e818bb6ff30370f8da86f5c092f27%26chksm%3Dce2d6f33f95ae625c18f68530a3d0af1f153d9719f771b8b3b8d55ff21433140b0e49e613d40%26scene%3D21%23wechat_redirect)

[全基因组关联分析（GWAS）扫不出信号怎么办（文献解读）](https://link.zhihu.com/?target=http%3A//mp.weixin.qq.com/s%3F__biz%3DMzg2MDA2MDQzMQ%3D%3D%26mid%3D2247483694%26idx%3D1%26sn%3Db7e1d4964b9edf9b8f8968d17d8d57fa%26chksm%3Dce2d6f4ef95ae65876b299ae1da014c490421598d41bb0830d0aa42b26db5fd7c523035a70a4%26scene%3D21%23wechat_redirect)

## 补充

GWAS其他教程：

www.transplantdb.eu/sites/transplantdb.eu/files/HandsOnTutorialtoGWAS\_Seren-030715.pdf

[https://doc.goldenhelix.com/SVS/tutorials/snp\_gwas/index.html](https://link.zhihu.com/?target=https%3A//doc.goldenhelix.com/SVS/tutorials/snp_gwas/index.html)

[http://ccbb.jnu.ac.in/IUBDDJan2015/workshop\_files/GWAS](https://link.zhihu.com/?target=http%3A//ccbb.jnu.ac.in/IUBDDJan2015/workshop_files/GWAS) Tutorial.pdf

[https://www.r-project.org/conferences/useR-2009/slides/Zhao+Tan.pdf](https://link.zhihu.com/?target=https%3A//www.r-project.org/conferences/useR-2009/slides/Zhao%2BTan.pdf)

users.du.se/~lrn/NOVAComputerExercises/NOVA\_GenABEL\_tutorial.pdf

[http://gsea4gwas-v2.psych.ac.cn/docs/tutorial.jsp](https://link.zhihu.com/?target=http%3A//gsea4gwas-v2.psych.ac.cn/docs/tutorial.jsp)

www.montefiore.ulg.ac.be/~kvansteen/GeneticEpi-UA2/Class5/Introduction to GenABEL.pdf

看看文献，加深对GWAS的理解：

[A tutorial on conducting genome‐wide association studies: Quality control and statistical analysis](https://link.zhihu.com/?target=https%3A//onlinelibrary.wiley.com/doi/full/10.1002/mpr.1608)

[Genome-wide association studies and beyond](https://link.zhihu.com/?target=https%3A//www.ncbi.nlm.nih.gov/pubmed/20235850)

[Genome-wide association studies](https://link.zhihu.com/?target=https%3A//www.ncbi.nlm.nih.gov/pubmed/23300413)

