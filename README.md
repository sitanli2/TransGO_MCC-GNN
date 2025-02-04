# TransGO Framework
The experimental source code of TransGO mentioned in the article "**Predicting Protein Function by Integrating Multi-Level Protein Sequence and Structural Features on a document based Database"** , along with the data scraping and splitting scripts involved.
![TransGO_refine3(英文版)](G:\360MoveData\Users\pc\Desktop\TransGo\Pics\TransGO模型图解构\整体结构\TransGO_refine3(英文版).png)
## 1. Dependencies
## 2. Dataset Preprocessing
### 2.1 ProGO Database
我们从一个文档型蛋白质关联数据库ProGO来获取我们所需的实验数据集。 具体的数据库搭建以及获取方式见 (url:)
- ProGO_scrapy/ 文件夹下的Uniprot_Scrapy, PDB_Scrapy, Chembl_Scrapy被用于从Uniprot, PDB, AlphaFold, Chembl 数据库的开放数据接口中获取蛋白质关联数据。并对接口返回的Json和XML文档进行数据清洗。
- ProGO_scrapy/MySQL_chembl 用于将ChEMBL32 的MySQL类型数据库转换为MonGODB可以存储的JSON字典。
### 2.2 Experimental Dataset PreProcessing
To predict protein functions by integrating sequence and structural data, we developed a comprehensive workflow based on ProGO database. Users can selected protein entries from ProGO (47,789 entries) based on specific criteria.
- Preprocessing/PDBID_mapping_UPID.py 用来根据指定蛋白质链的PDBID映射其在Uniprot数据库中的ID
- Preprocessing/CollectFromProGO.py 用来按照用户提供的筛选标准从ProGO中筛选出对应蛋白质的：完整序列，结构文件，功能标签
- Preprocessing/CollectFromDisk.py 用来精炼直接下载至本地的DeepFRI以及 PDB-B数据集
最后生成的数据集将被存储至UltraDataset, 其中的filteredIDlist.txt 为从ProGO中筛选出来的条目ID
如"1UMK$A=27-301$P00387" 其中：
- $P00387 代表该蛋白质条目的UPID
- 1UMK 为该蛋白质条目的PDBID
- $A 代表选择A链代表的空间结构映射
- =27-301代表A链在完整序列中的起始和结束位置
### 2.3 **Multi-level** Modeling of Protein Structures
We computes residue contact maps with different connection densities through Atoms(C) distance calculations and a unique hidden edge expansion method, calculating contact relationships including both explicit and implicit interactions.
- biotoolbox/ContactMapGen.py/contact_mapping_PDB()  用来从由实验测定的PDB文件中，按照不同接触阈值提取对应蛋白质链的contact map
- biotoolbox/ContactMapGen.py/contact_mapping_AlphaFold() 用来从由AlphaFold预测的结构文件中，按照不同接触阈值提取对应蛋白质链的contact map
- biotoolbox/ContactMapGen.py/ContactMap_densities_Switching() 用来寻找强相关性的隐藏边，或者是弱相关性的冗余边。
- biotoolbox/ContactMapGen.py/insert_hidden_edge() 用来给对应的接触图进行隐藏边扩充
- biotoolbox/ContactMap_visualization.py 用于蛋白质接触图可视化
- biotoolbox/AlphaFoldfile_Dataset_Ultimate.zip and PDBfile_Dataset_Ultimate.zip 为初始pdb结构文件
- biotoolbox/Contactmapfile_Dataset_Ultimate.zip 为提取出的有不同接触密度的蛋白质接触图
## 3. Train and validation
We conduct a series of controlled and ablation experiments to evaluate the performance of the TransGO model and to investigate methods for reducing training costs.
- TransGO_refine/DeepFRI_main.py and ESM_main.py: 基于MaW et al and Gligorijević V et al等人的工作复现DeepFRI模型
- ESM_extract and TransGO_refine/utils/ESM_bilstm_Dataset():  combined two sophisticated protein language models, the model introduced by Belpler et.al(**Url**:) based on Bi-LSTM, and the Evolutionary Scale Modeling (ESM) based on Transformer(**Url**:), to create a hybrid feature extractor.
- TransGO_refine/model.py 集成了多种不同的图神经网络架构，其中包括 MCC-GNN :A multi-channel Graph Neural Networks that integrates sequence representations and contact maps with varying contact densities. The network features three configurations—single, dual, and triple channels, each offering unique ways to capture and combine different aspects of protein data.
- TransGO_refine/TransGO_main.py 对TransGO模型以及其它对照模型进行训练与验证 （在训练的过程中更换双通道MCC-GNN不同通道的残基接触图）
