# 选取训练集中最常见的10个 最常见的GO 条目，并且计算了它们的占比.然后对测试集中的每个蛋白质条目都做同样的预测.
import numpy as np
from utils import get_results
from utils import get_label
from sklearn.metrics import precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt
import copy


def get_GOlabel(func = 'mf'):
    if func == 'mf':
        mf_dict = {}
        mf_func = []
        with open('../UltraDataset/Database1+/Go_label.txt', 'r') as f:
            for line in f:
                if '>' in line:
                    pdb_chain_uid = line[1:].strip('\n')
                elif 'mf' in line:
                    if 'mf:' == line.strip('\n'):
                        continue
                    else:
                        mf_func_list = line[3:].strip().split('\t')
                        mf_dict[pdb_chain_uid] = mf_func_list
        with open('../UltraDataset/Database1+/mf/mf_label.txt', 'r') as f:
            for line in f:
                line = line.strip('\n')
                mf_func.append(line)

        return mf_dict,len(mf_func),mf_func

    elif func == 'bp':
        bp_dict = {}
        bp_func = []
        with open('../UltraDataset/Database1+/Go_label.txt', 'r') as f:
            for line in f:
                if '>' in line:
                    pdb_chain_uid = line[1:].strip('\n')
                elif 'bp' in line:
                    if 'bp:' == line.strip('\n'):
                        continue
                    else:
                        bp_func_list = line[3:].strip().split('\t')
                        bp_dict[pdb_chain_uid] = bp_func_list
        with open('../UltraDataset/Database1+/bp/bp_label.txt', 'r') as f:
            for line in f:
                line = line.strip('\n')
                bp_func.append(line)

        return bp_dict, len(bp_func), bp_func
    elif func == 'cc':
        cc_dict = {}
        cc_func = []
        with open('../UltraDataset/Database1+/Go_label.txt', 'r') as f:
            for line in f:
                if '>' in line:
                    pdb_chain_uid = line[1:].strip('\n')
                elif 'cc' in line:
                    if 'cc:' == line.strip('\n'):
                        continue
                    else:
                        cc_func_list = line[3:].strip().split('\t')
                        cc_dict[pdb_chain_uid] = cc_func_list
        with open('../UltraDataset/Database1+/cc/cc_label.txt', 'r') as f:
            for line in f:
                line = line.strip('\n')
                cc_func.append(line)

        return cc_dict, len(cc_func), cc_func

def GetTop10_GO(function = 'mf'):
    """得到前十个最频繁出现的GO标签"""

    label_dict,funclist_len,funclist = get_GOlabel(func= function) #得到 mf 的label_dict
    # print(label_dict,'\n',funclist_len,'\n',funclist)
    GOfrequency_dict = {}
    for GO in funclist:
        GOfrequency_dict[GO] = 0 #初始化frequency 字典，全部为零
    # print(label_dict,'\n',len(label_dict))
    # print(GOfrequency_dict,'\n',len(GOfrequency_dict))

    # count = 0
    # for key, value in label_dict.items():
    #     count += len(value)
    # print(count)

    for key1, value1 in label_dict.items():
        for GOitem in value1:
            GOfrequency_dict[GOitem] = GOfrequency_dict[GOitem] + 1
    print(GOfrequency_dict,'\n',len(GOfrequency_dict))

    """生成代表所有GOterms 在数据集中出现次数（频率）的frequency字典后
    按照出现频数给所有GO从大到小排序，生成新的排序后的字典，而且，取排序后字典的前10个GO id生成一个列表,作为naive的分配值
    """
    # 按照值（value）从大到小排序frequency字典，并转换为列表
    sorted_items = sorted(GOfrequency_dict.items(), key=lambda x: x[1], reverse=True)

    # 创建排序后的frequency新字典
    sorted_dict = dict(sorted_items)
    print(sorted_dict, '\n', len(sorted_dict))
    # 提取排序后字典的前10的 GO label
    top_10_keys = [key for key, value in sorted_items[:10]]
    print(f'top_10_keys of the training set on {function} label:',top_10_keys)
    return top_10_keys

def GetpreDict_TrueDict(GO_label_list,functionset = 'mf'):
    trueDict, labelnumber = get_label(func = functionset)  # 真实标签
    preDict = {}
    for key, values in trueDict.items():
        # 创建一个维度为10，所有值为1的numpy数组
        preDict[key] = np.ones(10)
    # print(preDict)

    New_trueDict = copy.deepcopy(preDict)  # 制作新的只考虑有没有这TOP10个GO标签的dict
    label_dict, funclist_len, funclist = get_GOlabel(func = functionset)  # 得到 mf 的label_dict
    # print(label_dict)

    for key1, values1 in New_trueDict.items():
        i = 0  # 10维numpy数组的定位指针
        true_label_list = label_dict[key1]
        for item in GO_label_list:
            if (item in true_label_list):  # 此时说明确实存在这个GO标签,不需要置零
                i += 1
                continue
            else:  # 没这个GO terms,需要置零
                values1[i] = 0
                i += 1
    # print(New_trueDict)
    # print(preDict)
    return New_trueDict,preDict

def AUPR_test(true_dict,pred_dict):
    # 示例数据：真实标签和预测标签（根据你的描述，这里用的是假设的数据）
    # 假设 true_dict 和 pred_dict 的结构如下：
    # true_dict = {
    #     '5ZJR-A_P09087': np.array([0., 0., 0., 1., 0., 1., 0., 0., 1., 0.]),
    #     '5ZJR-A_P09088': np.array([1., 0., 0., 1., 0., 0., 1., 0., 1., 0.]),
    #     # 更多样本...
    # }
    #
    # pred_dict = {
    #     '5ZJR-A_P09087': np.array([0.1, 0.2, 0.3, 0.9, 0.4, 0.8, 0.2, 0.3, 0.7, 0.5]),
    #     '5ZJR-A_P09088': np.array([0.9, 0.1, 0.3, 0.8, 0.4, 0.2, 0.7, 0.1, 0.6, 0.9]),
    #     # 更多样本...
    # }

    # 提取所有样本的真实标签和预测标签
    # 获取所有的样本ID
    sample_ids = list(true_dict.keys())

    # 初始化真实标签矩阵和预测标签矩阵
    y_true = []
    y_pred = []

    for sample_id in sample_ids:
        y_true.append(true_dict[sample_id])  # 真实标签
        y_pred.append(pred_dict[sample_id])  # 预测标签

    # 转换为 NumPy 数组，方便后续计算
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    # 计算每个标签的 Precision-Recall 曲线和 AUPR
    precision = {}
    recall = {}
    f_max = {}
    average_precision = {}

    for i in range(y_true.shape[1]):  # 遍历每个GO标签
        precision[i], recall[i], _ = precision_recall_curve(y_true[:, i], y_pred[:, i])
        average_precision[i] = average_precision_score(y_true[:, i], y_pred[:, i])
        # 计算 F1 分数（不同阈值下的 F1）
        f1_scores = 2 * (precision[i] * recall[i]) / (precision[i] + recall[i])
        # F-max 是 F1 分数的最大值
        f_max[i] = np.max(f1_scores)

    # 绘制 PR 曲线
    plt.figure(figsize=(10, 6))
    for i in range(y_true.shape[1]):
        plt.plot(recall[i], precision[i], label=f'GO Label {i + 1} (AUPR={average_precision[i]:.2f})')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve for each GO Label')
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()

    # # 打印每个标签的 AUPR
    # print("AUPR for each label:")
    # for i in range(y_true.shape[1]):
    #     print(f"GO Label {i + 1}: {average_precision[i]:.2f}")

    # 打印每个标签的 AUPR 和 F-max
    print("AUPR and F-max for each label:")
    for i in range(y_true.shape[1]):
        print(f"GO Label {i + 1}: AUPR={average_precision[i]:.2f}, F-max={f_max[i]:.2f}")

    # 计算总体 AUPR（可以将所有标签的 AUPR 值汇总，计算加权平均）
    mean_aupr = np.mean(list(average_precision.values()))
    mean_f_max = np.mean(list(f_max.values()))
    print(f"Mean AUPR: {mean_aupr:.2f}")
    print(f"Mean F-max: {mean_f_max:.2f}")




if __name__ == '__main__':
    """得到前十个最频繁出现的GO标签"""
    # GO_label_list = GetTop10_GO(function='bp') # 得到前十个最频繁出现的GO标签,如下：
    GO_label_list_mf = ['GO:0097367', 'GO:0032553', 'GO:0017076', 'GO:0032555', 'GO:0035639', 'GO:0140096', 'GO:0030554', 'GO:0032559', 'GO:0005524', 'GO:0003677']
    GO_label_location_mf = [245,449,462,117,302,225,293,155,14,326] #最常出现的前十个mf标签在共486个标签mf_label.txt中的位置
    GO_label_list_bp = ['GO:0051716', 'GO:0071840', 'GO:0019222', 'GO:0016043', 'GO:0016070', 'GO:0019438', 'GO:0060255', 'GO:0051179', 'GO:0031323', 'GO:0080090']
    GO_label_list_cc = ['GO:0005634', 'GO:0005829', 'GO:0031974', 'GO:0070013', 'GO:0043233', 'GO:0071944', 'GO:0043228', 'GO:0043232', 'GO:0031981', 'GO:0005886']


    """分配标签给测试集(考虑全486个标签)并计算AUPR值"""
    # trueDict, labelnumber = get_label(func='mf')  # 真实标签
    # # print(trueDict)
    # preDict = copy.deepcopy(trueDict)
    # for key,values in preDict.items():
    #     # 将现有numpy数组的所有元素重置为0
    #     values.fill(0)
    #
    #     for loc in GO_label_location:
    #         # 将第loc维（索引位置为loc-1）置为1
    #         values[loc-1] = 1
    #
    # # print(trueDict,'\n',preDict)
    # AUPR_test(trueDict,preDict)

    """分配标签给测试集(仅考虑前10个标签)并计算AUPR值"""
    New_trueDict,preDict = GetpreDict_TrueDict(GO_label_list_bp,functionset='bp')
    # print(New_trueDict)
    # print(preDict)
    AUPR_test(New_trueDict, preDict)










