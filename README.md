# kmldm
不同环境条件下的关联网络分析
## Split process
首先按照环境因素对样本进行分割, 每次将数据进行2-GMM, 递归进行直到cluster中的样本数少于某个阈值
## Merge process
对split产生的叶子节点进行关联网络分析, 根据估计出的参数进行节点合并
