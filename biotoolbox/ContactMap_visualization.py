import numpy as np
import pickle as pk
import matplotlib.pyplot as plt

def test_on_NPYfile(npy_file):
    test_array = np.load(npy_file)
    plt.figure(figsize=(10, 5))  # 这行代码创建了一个新的图形窗口，并指定了该图形的大小为宽度10英寸、高度5英寸。
    plt.subplot(1, 2, 1)
    plt.imshow(test_array, cmap='gray')
    plt.title('grey pic')

    plt.subplot(1, 2, 2)  # 选择第二个子图进行绘制
    plt.imshow(test_array, cmap='viridis')  # viridis 色彩映射，它是一种颜色渐变映射。
    plt.title('threshold=')
    plt.show()


def test_on_NPZfile():
    data = np.load('G:/fkfk/Protein_Predict/DeepFRI-master/examples/pdb_cmaps/1S3P-A.npz')
    '''打印npz文件中的键（keys）确定其中包含数组的名称'''
    print("keys in the compare_array 1S3P-A.npz file:",data.files,'\n') #['C_alpha', 'C_beta', 'seqres']
    #接下来提取需要做可视化的数组，如果数组是图像数据，可以使用matplotlib进行可视化
    npz_arry1 = data['C_alpha']
    npz_arry2 = data['C_beta']
    npz_arry3 = data['seqres']

    print('npz_arry1即 C_alpha：',npz_arry1,'\n')
    print('npz_arry2即 C_beta：',npz_arry2,'\n')
    print('npz_arry3即 seqres：',npz_arry3,'\n')

    plt.figure(figsize=(10,5)) #这行代码创建了一个新的图形窗口，并指定了该图形的大小为宽度10英寸、高度5英寸。
    plt.subplot(1, 3, 1)
    '''
    plt.subplot(1, 2, 1) 
    这行代码创建了一个包含 1 行 3 列的子图网格，并选择了第一个子图来绘制。
    参数 (1, 3, 1) 指定了子图网格的布局，其中 (1, 3) 表示总共有1行3列子图，而 ，1） 表示当前选中的是第一个子图。
    '''
    plt.imshow(npz_arry1, cmap='gray')
    '''
    plt.imshow(npz_arry1, cmap='gray')
    这行代码使用 imshow() 函数将名为 npz_arry1 的 NumPy 数组绘制成图像 其中:
    - npz_arry1 是一个二维数组，通常是表示图像的像素值矩阵。
    - imshow() 函数用于将二维数组中的数值映射为颜色，然后以图像的形式显示出来。
    - cmap='gray' 指定了使用灰度色彩映射（colormap），即将二维数组中的值映射到灰度颜色空间中，显示出黑白图像。这意味着数值较低的点将显示为较暗的灰色，而数值较高的点将显示为较亮的灰色。
    '''
    plt.title('npz_arry1: C_alpha')

    plt.subplot(1, 3, 2) # 选择第二个子图进行绘制
    plt.imshow(npz_arry2, cmap='viridis') # viridis 色彩映射，它是一种颜色渐变映射。
    plt.title('npz_arry2: C_beta')

    # plt.subplot(1, 3, 3)
    # plt.imshow(npz_arry3, cmap='jet') # cmap='jet' 指定了使用 jet 色彩映射，它是一种带有彩虹色彩的映射，适用于显示多种颜色。
    # plt.title('5XYI-g.npy')


    plt.show()

def drawcontactmap(file_source_list):
    plt.figure(figsize=(10, 6))  # 这行代码创建了一个新的图形窗口，并指定了该图形的大小为宽度8英寸、高度6英寸。
    i=1
    thre = 6
    for contact_martrix_file_source in file_source_list:
        test_array = np.load(contact_martrix_file_source)
        # plt.subplot(1, 8, i) # 定义8个子图，这是第i个
        # plt.imshow(test_array, cmap='gray')
        # plt.title('threshold='+thre+"gray")
        # i += 1

        # 创建一个2行2列的子图
        # 第一个子图 (i=1,位于第 1 行第 1 列)
        plt.subplot(2, 2, i)  # 再画一个彩色的子图
        plt.imshow(test_array, cmap='viridis')  # viridis 色彩映射，它是一种颜色渐变映射。
        plt.title('threshold='+str(thre))
        i += 1
        thre += 2
        # 调整子图之间的间距
    plt.tight_layout()
    plt.show()


def drawcontactmap_weiyi(file_source,title_thre):
    plt.figure(figsize=(8, 6))
    test_array = np.load(file_source)
    plt.imshow(test_array, cmap='viridis')  # viridis 色彩映射，它是一种颜色渐变映射。
    # plt.title('threshold=' + str(title_thre))
    plt.tight_layout()
    plt.show()

def evaluation_visualization_upgrade():
    """
    使用python，画一个2*3的图，

    前三张图为折线图，分别代表method1,mothod2,method3的F-max评分折线，其中method1折线上数据点为三角形，method2为圆形，method3为正方形。
    横轴为epoch: [1, 5, 10, 20, 30, 40, 50]，纵轴为F-MAX: [0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8].
    其中method1的在epoch: [1, 5, 10, 20, 30, 40, 50] 的对应F-max列表为：[0.212, 0.462, 0.576, 0.657, 0.678, 0.684, 0.711]..method2,method3也类似

    后三张图为柱状图，分别代表method1,mothod2,method3的AUPR评分,其中method1的柱状图颜色为蓝色，method2为黄色，method3为红色。
    横轴epoch:[10, 20, 30, 40, 50],纵轴AUPR：[0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]，
    其中method1的在epoch: [10, 20, 30, 40, 50] 的对应AUPR列表为：[0.584, 0.685, 0.706, 0.710, 0.715]....method2,method3也类似

    :return:
    """
    # 准备数据
    epochs = [1, 5, 10, 20, 30, 40, 50]
    epochs2 = [10, 20, 30, 40, 50]
    aupr_scores_method1_MF = [0.168, 0.453, 0.584, 0.685, 0.706, 0.710, 0.715]  # DeepFRI的m-AUPR分数
    f1_scores_method1_MF = [0.212, 0.462, 0.576, 0.657, 0.678, 0.684, 0.711]  # DeepFRI的F1分数
    aupr_scores_method2_MF = [0.199, 0.467, 0.616, 0.710, 0.723, 0.729, 0.734]  # TransGO的AUPR分数
    f1_scores_method2_MF = [0.236, 0.480, 0.603, 0.690, 0.697, 0.707, 0.714]  # TransGO的F1分数
    aupr_scores_method3_MF = []  # DeepFRI(ESM2)的m-AUPR分数
    f1_scores_method3_MF = []  # DeepFRI(ESM2)的F1分数

    aupr_scores_method1_BP = [0.10, 0.167, 0.267, 0.353, 0.386, 0.401, 0.407]  # DeepFRI的AUPR分数
    f1_scores_method1_BP = [0.172, 0.243, 0.318, 0.385, 0.411, 0.422, 0.435]  # DeepFRI的F1分数
    aupr_scores_method2_BP = [0.123, 0.217, 0.292, 0.364, 0.399, 0.405, 0.434]  # TransGO的AUPR分数
    f1_scores_method2_BP = [0.191, 0.282, 0.334, 0.389, 0.423, 0.438, 0.454]  # TransGO的F1分数
    aupr_scores_method3_BP = []  # DeepFRI(ESM2)的m-AUPR分数
    f1_scores_method3_BP = []  # DeepFRI(ESM2)的F1分数

    aupr_scores_method1_CC = [0.120, 0.340, 0.445, 0.478, 0.495, 0.506, 0.509]  # DeepFRI的AUPR分数
    f1_scores_method1_CC = [0.155, 0.387, 0.471, 0.496, 0.512, 0.513, 0.521]  # DeepFRI的F1分数
    aupr_scores_method2_CC = [0.136, 0.358, 0.469, 0.485, 0.501, 0.510, 0.511]  # TransGO的AUPR分数
    f1_scores_method2_CC = [0.177, 0.390, 0.482, 0.511, 0.515, 0.519, 0.527]  # TransGO的F1分数
    aupr_scores_method3_CC = []  # DeepFRI(ESM2)的m-AUPR分数
    f1_scores_method3_CC = []  # DeepFRI(ESM2)的F1分数

    method1_aupr = [0.2, 0.3, 0.4, 0.5, 0.6]  # 示例数据，替换为您的实际数据
    method2_aupr = [0.25, 0.35, 0.45, 0.55, 0.65]  # 示例数据，替换为您的实际数据
    method3_aupr = [0.3, 0.4, 0.5, 0.6, 0.7]  # 示例数据，替换为您的实际数据

    # 创建子图
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # 图表标题
    titles = ["MF", "BP", "CC"]
    for i, ax in enumerate(axes):
        for j, ax2 in enumerate(ax):
            if (i == 0): # 第一行的三个子图 (三个折线图)
                if(titles[j] == 'MF'):
                    ax2.plot(epochs, aupr_scores_method1_MF, 'g--^', label="DeepFRI m-AUPR")  # 绿色虚线+三角形 (method1, AUPR)
                    ax2.plot(epochs, f1_scores_method1_MF, 'g--s', label="DeepFRI F1-score")  # 绿色虚线+方形 (method1, F1-score)
                    ax2.plot(epochs, aupr_scores_method2_MF, 'b-^', label="TransGO m-AUPR")  # 蓝色实线+三角形 (method2, AUPR)
                    ax2.plot(epochs, f1_scores_method2_MF, 'b-s', label="TransGO F1-score")  # 蓝色实线+方形 (method2, F1-score)
                elif (titles[j] == 'BP'):
                    ax2.plot(epochs, aupr_scores_method1_BP, 'g--^', label="DeepFRI m-AUPR")  # 绿色虚线+三角形 (method1, AUPR)
                    ax2.plot(epochs, f1_scores_method1_BP, 'g--s',label="DeepFRI F1-score")  # 绿色虚线+方形 (method1, F1-score)
                    ax2.plot(epochs, aupr_scores_method2_BP, 'b-^', label="TransGO m-AUPR")  # 蓝色实线+三角形 (method2, AUPR)
                    ax2.plot(epochs, f1_scores_method2_BP, 'b-s',label="TransGO F1-score")  # 蓝色实线+方形 (method2, F1-score)
                elif (titles[j] == 'CC'):
                    ax2.plot(epochs, aupr_scores_method1_CC, 'g--^', label="DeepFRI m-AUPR")  # 绿色虚线+三角形 (method1, AUPR)
                    ax2.plot(epochs, f1_scores_method1_CC, 'g--s',label="DeepFRI F1-score")  # 绿色虚线+方形 (method1, F1-score)
                    ax2.plot(epochs, aupr_scores_method2_CC, 'b-^', label="TransGO m-AUPR")  # 蓝色实线+三角形 (method2, AUPR)
                    ax2.plot(epochs, f1_scores_method2_CC, 'b-s',label="TransGO F1-score")  # 蓝色实线+方形 (method2, F1-score)
                ax2.set_title(titles[j])  # 设置子图标题
                ax2.set_xlabel("Epoch")  # 设置x轴标签
                ax2.set_ylabel("Evaluation Score")  # 设置y轴标签
                ax2.set_xticks(epochs)  # 设置x轴刻度
                ax2.set_yticks([0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8])  # 设置y轴刻度
                ax2.legend(loc="lower right")  # 添加图例

            elif (i == 1):# 第二行的三个子图 (三个柱状图)
                # 设置柱状图宽度
                width = 0.2

                # 设置x轴位置
                x = np.arange(len(epochs))


                ax2.bar(x - width, method1_aupr, width, label='Method 1', color='skyblue')
                ax2.bar(x, method2_aupr, width, label='Method 2', color='lightgreen')
                ax2.bar(x + width, method3_aupr, width, label='Method 3', color='lightcoral')

                # 设置子图标题和标签
                ax2.set_title(titles[j])
                ax2.set_xlabel('Epoch')
                ax2.set_ylabel('AUPR')

                # 设置x轴刻度和标签
                ax2.set_xticks(x)
                ax2.set_xticklabels(epochs)

                # 设置y轴刻度
                ax2.set_yticks(np.arange(0, 0.81, 0.05))

                # 添加图例
                ax2.legend()

                # 添加数值标签
                for j, v in enumerate(method1_aupr):
                    ax.text(j - width, v + 0.01, str(v), ha='center', va='bottom', fontsize=8)
                for j, v in enumerate(method2_aupr):
                    ax.text(j, v + 0.01, str(v), ha='center', va='bottom', fontsize=8)
                for j, v in enumerate(method3_aupr):
                    ax.text(j + width, v + 0.01, str(v), ha='center', va='bottom', fontsize=8)

    # 显示图形
    plt.tight_layout()
    # plt.savefig("G:/360MoveData/Users/pc/Desktop/evaluation_scores.png", dpi=300)  # 保存图片，设置分辨率为300dpi
    plt.show()


aupr_scores_method1_MF = [0.168, 0.453, 0.584, 0.685, 0.706, 0.710, 0.715]  # DeepFRI的m-AUPR分数
aupr_scores_method2_MF = [0.199, 0.467, 0.616, 0.710, 0.723, 0.729, 0.734]  # TransGO的AUPR分数
aupr_scores_method3_MF = []  # DeepFRI(ESM2)的m-AUPR分数

aupr_scores_method1_BP = [0.10, 0.167, 0.267, 0.353, 0.386, 0.401, 0.407]  # DeepFRI的AUPR分数
aupr_scores_method2_BP = [0.123, 0.217, 0.292, 0.364, 0.399, 0.405, 0.434]  # TransGO的AUPR分数
aupr_scores_method3_BP = []  # DeepFRI(ESM2)的m-AUPR分数

aupr_scores_method1_CC = [0.120, 0.340, 0.445, 0.478, 0.495, 0.506, 0.509]  # DeepFRI的AUPR分数
aupr_scores_method2_CC = [0.136, 0.358, 0.469, 0.485, 0.501, 0.510, 0.511]  # TransGO的AUPR分数
aupr_scores_method3_CC = []  # DeepFRI(ESM2)的m-AUPR分数


def evaluation_visualization():
    # 准备数据
    epochs = [1, 5, 10, 20, 30, 40, 50]
    # epochs2 = [10, 20, 30, 40, 50]

    f1_scores_method1_MF = [0.212, 0.462, 0.575, 0.656, 0.677, 0.680, 0.684]  # DeepFRI的F1分数
    f1_scores_method2_MF = [0.236, 0.480, 0.603, 0.690, 0.697, 0.713, 0.718]  # TransGO的F1分数
    f1_scores_method3_MF = [0.221, 0.413, 0.486, 0.533, 0.542, 0.541, 0.549]  # DeepFRI(ESM2)的F1分数


    f1_scores_method1_BP = [0.172, 0.237, 0.300, 0.367, 0.403, 0.420, 0.421]  # DeepFRI的F1分数
    f1_scores_method2_BP = [0.140, 0.247, 0.324, 0.389, 0.413, 0.433, 0.437]  # TransGO的F1分数
    f1_scores_method3_BP = [0.190, 0.241, 0.291, 0.308, 0.324, 0.326, 0.330]  # DeepFRI(ESM2)的F1分数


    f1_scores_method1_CC = [0.154, 0.387, 0.471, 0.489, 0.502, 0.499, 0.503]  # DeepFRI的F1分数
    f1_scores_method2_CC = [0.199, 0.412, 0.482, 0.511, 0.520, 0.518, 0.523]  # TransGO的F1分数
    f1_scores_method3_CC = [0.217, 0.402, 0.438, 0.440, 0.441, 0.439, 0.446]  # DeepFRI(ESM2)的F1分数

    # 创建子图
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 图表标题
    titles = ["MF", "BP", "CC"]

    # 绘制每个子图
    for i, ax in enumerate(axes):
        if (titles[i] == 'MF'):
            # ax.plot(epochs, aupr_scores_method1_MF, 'g--s', label="DeepFRI m-AUPR")  # 绿色虚线+方形 (method1, AUPR)
            ax.plot(epochs, f1_scores_method1_MF, 'g-^', label="DeepFRI")  # 绿色实线+三角形 (method1, F1-score)
            # ax.plot(epochs, aupr_scores_method2_MF, 'b-^', label="TransGO m-AUPR")  # 蓝色实线+三角形 (method2, AUPR)
            ax.plot(epochs, f1_scores_method2_MF, 'r-^', label="TransGO")  # 红色实线+三角形 (method2, F1-score)
            ax.plot(epochs, f1_scores_method3_MF, 'b-^', label="DeepFRI(ESM2)")  # 蓝色实线+三角形 (method3, F1-score)
        elif (titles[i] == 'BP'):
            # ax.plot(epochs, aupr_scores_method1_BP, 'g--^', label="DeepFRI m-AUPR")  # 绿色虚线+三角形 (method1, AUPR)
            ax.plot(epochs, f1_scores_method1_BP, 'g-^', label="DeepFRI")  # 绿色虚线+方形 (method1, F1-score)
            # ax.plot(epochs, aupr_scores_method2_BP, 'b-^', label="TransGO m-AUPR")  # 蓝色实线+三角形 (method2, AUPR)
            ax.plot(epochs, f1_scores_method2_BP, 'r-^', label="TransGO")  # 红色实线+三角形 (method2, F1-score)
            ax.plot(epochs, f1_scores_method3_BP, 'b-^', label="DeepFRI(ESM2)")  # 蓝色实线+方形 (method3, F1-score)
        elif (titles[i] == 'CC'):
            # ax.plot(epochs, aupr_scores_method1_CC, 'g--^', label="DeepFRI m-AUPR")  # 绿色虚线+三角形 (method1, AUPR)
            ax.plot(epochs, f1_scores_method1_CC, 'g-^', label="DeepFRI")  # 绿色虚线+方形 (method1, F1-score)
            # ax.plot(epochs, aupr_scores_method2_CC, 'b-^', label="TransGO m-AUPR")  # 蓝色实线+三角形 (method2, AUPR)
            ax.plot(epochs, f1_scores_method2_CC, 'r-^', label="TransGO")  # 红色实线+三角形 (method2, F1-score)
            ax.plot(epochs, f1_scores_method3_CC, 'b-^', label="DeepFRI(ESM2)")  # 蓝色实线+方形 (method3, F1-score)

        ax.set_title(titles[i])  # 设置子图标题
        ax.set_xlabel("Epoch")  # 设置x轴标签
        ax.set_ylabel("F-MAX Score")  # 设置y轴标签
        ax.set_xticks(epochs)  # 设置x轴刻度
        ax.set_yticks([0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8])  # 设置y轴刻度
        ax.legend(loc="lower right")  # 添加图例

    # 显示图形
    plt.tight_layout()
    plt.savefig("G:/360MoveData/Users/pc/Desktop/evaluation_scores.png", dpi=300)  # 保存图片，设置分辨率为300dpi
    plt.show()

def zhuzhuangtu():
    # 数据
    epochs = [10, 20, 30, 40, 50]

    method1_aupr = [0.2, 0.3, 0.4, 0.5, 0.6]  # 示例数据，请替换为您的实际数据
    method2_aupr = [0.25, 0.35, 0.45, 0.55, 0.65]  # 示例数据，请替换为您的实际数据
    method3_aupr = [0.3, 0.4, 0.5, 0.6, 0.7]  # 示例数据，请替换为您的实际数据

    aupr_scores_method1_MF = [0.583, 0.685, 0.706, 0.709, 0.712]  # DeepFRI的m-AUPR分数
    aupr_scores_method2_MF = [0.616, 0.710, 0.723, 0.729, 0.734]  # TransGO的AUPR分数
    aupr_scores_method3_MF = [0.486, 0.554, 0.542, 0.548, 0.562]  # DeepFRI(ESM2)的m-AUPR分数

    aupr_scores_method1_BP = [0.235, 0.323, 0.373, 0.390, 0.391]  # DeepFRI的AUPR分数
    aupr_scores_method2_BP = [0.248, 0.348, 0.386, 0.399, 0.406]  # TransGO的AUPR分数
    aupr_scores_method3_BP = [0.222, 0.261, 0.275, 0.282, 0.284]  # DeepFRI(ESM2)的m-AUPR分数

    aupr_scores_method1_CC = [0.398, 0.462, 0.480, 0.481, 0.483]  # DeepFRI的AUPR分数
    aupr_scores_method2_CC = [0.440, 0.505, 0.505, 0.504, 0.509]  # TransGO的AUPR分数
    aupr_scores_method3_CC = [0.395, 0.389, 0.404, 0.419, 0.422]  # DeepFRI(ESM2)的m-AUPR分数

    # 子图标题
    titles = ['MF', 'BP', 'CC']

    # 创建子图
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # 1行3列的子图

    # 设置柱状图宽度
    width = 0.28

    # 设置x轴位置
    x = np.arange(len(epochs))

    # 循环绘制子图
    for i, ax in enumerate(axes):
        # 绘制柱状图
        if (titles[i] == 'MF'):
            method1_aupr = aupr_scores_method1_MF
            method2_aupr = aupr_scores_method2_MF
            method3_aupr = aupr_scores_method3_MF
        elif (titles[i] == 'BP'):
            method1_aupr = aupr_scores_method1_BP
            method2_aupr = aupr_scores_method2_BP
            method3_aupr = aupr_scores_method3_BP
        elif (titles[i] == 'CC'):
            method1_aupr = aupr_scores_method1_CC
            method2_aupr = aupr_scores_method2_CC
            method3_aupr = aupr_scores_method3_CC

        ax.bar(x - width, method1_aupr, width, label='DeepFRI', color='lightgreen')
        ax.bar(x, method2_aupr, width, label='TransGO', color='lightcoral')
        ax.bar(x + width, method3_aupr, width, label='DeepFRI(ESM2)', color='skyblue')

        # 设置子图标题和标签
        ax.set_title(titles[i])
        ax.set_xlabel('Epoch')
        ax.set_ylabel('AUPR')

        # 设置x轴刻度和标签
        ax.set_xticks(x)
        ax.set_xticklabels(epochs)
        if (titles[i] == 'MF'):
            # 设置y轴刻度
            ax.set_yticks(np.arange(0, 1.01, 0.05))
        elif(titles[i] == 'BP' or titles[i] == 'CC'):
            ax.set_yticks(np.arange(0, 0.81, 0.05))

        # 添加图例
        ax.legend()

        # 添加数值标签（缩小字号）
        for j, v in enumerate(method1_aupr):
            ax.text(j - width, v + 0.01, f"{v:.3f}", ha='center', va='bottom', fontsize=6)  # 修改此处
        for j, v in enumerate(method2_aupr):
            ax.text(j, v + 0.01, f"{v:.3f}", ha='center', va='bottom', fontsize=6)  # 修改此处
        for j, v in enumerate(method3_aupr):
            ax.text(j + width, v + 0.01, f"{v:.3f}", ha='center', va='bottom', fontsize=6)  # 修改此处

    # 调整子图布局
    plt.tight_layout()
    plt.savefig("G:/360MoveData/Users/pc/Desktop/evaluation_scores_AUPR.png", dpi=300)  # 保存图片，设置分辨率为300dpi
    # 显示图形
    plt.show()


if __name__ == '__main__':
    # test_on_NPZfile()

    """接触图可视化"""
    # threshold_list = ['6','8','10','12']
    # file_source_list = []
    # for threshold in threshold_list:
    #     file_source_list.append('G:/fkfk/AlphaFold_refine/biotoolbox/Test_ContactMapFiles/1afs_'+threshold+'.0A-intact.npy') #6A,8A,10A,12A
    # # drawcontactmap(file_source_list)
    #
    # thre_test = "10"
    # drawcontactmap_weiyi(file_source='G:/fkfk/AlphaFold_refine/biotoolbox/Test_ContactMapFiles/1afs_'+thre_test+'.0A-intact.npy',title_thre=thre_test)

    """evaluation score visualization"""
    # evaluation_visualization_upgrade()
    evaluation_visualization()
    zhuzhuangtu()


'''
结论：
.npy 和 .npz 都是与 NumPy 库相关的文件扩展名，用于存储 NumPy 数组数据。它们之间的主要区别在于存储方式和文件结构：

.npy 文件：

.npy 文件是 NumPy 的二进制格式，用于存储单个 NumPy 数组。
它以二进制格式存储数组数据，并且可以保留数组的结构、形状、dtype 等信息。
由于是单个数组，.npy 文件适用于存储单个数组或加载单个数组。
.npz 文件：

.npz 文件也是 NumPy 的二进制格式，但是它可以存储多个 NumPy 数组。
.npz 文件实际上是一个压缩文件，其中包含了多个 .npy 格式的数组以及数组的名称。
这种格式适用于需要同时保存和加载多个 NumPy 数组的情况，例如在实验数据处理中，通常会产生多个相关联的数组。
综上所述，.npy 文件用于单个数组的存储，而 .npz 文件用于存储多个相关联的数组，并且能够将它们以压缩格式存储在一个文件中，方便管理和传输。

'''


'''
对npz 文件分析后可得出结论，每个npz文件包含三个部分:
keys in the compare_array 1S3P-A.npz file: ['C_alpha', 'C_beta', 'seqres'] 
其中C_alpha 和 C_beta 都为npy文件即numpy数组
而seqres 则为一字符串文件 代表当前如 1S3P-A.npz 对应的蛋白质序列
seqres： SMTDLLSAEDIKKAIGAFTAADSFDHKKFFQMVGLKKKSADDVKKVFHILDKDKDGFIDEDELGSILKGFSSDARDLSAKETKTLMAAGDKDGDGKIGVEEFSTLVAES 




'''

# def persistent_load(pers_id):
#     return None
#
# with open('G:/fkfk/Struct2Go-main/data_collect/amplify_samples/model/bp/30/pdb_amplifed_bp_GCN_512_gap-gmp_0.2_8.0.pkl','rb') as file1:
#     unpickler= pk.Unpickler(file1)
#     unpickler.persistent_load = persistent_load()
#
#     loaded_data = unpickler.load()
#
# with open('G:/fkfk/Struct2Go-main/data_collect/amplify_samples/model/bp/30/pdb_no_amplifed_bp_GCN_512_gap-gmp_0.2_8.0.pkl','rb') as file2:
#     unpickler = pk.Unpickler(file2)
#     unpickler.persistent_load = persistent_load()
#
#     loaded_data1 = unpickler.load()
#
# print("pdb_amplifed_bp_GCN_512_gap-gmp_0.2_8.0.pkl = ",loaded_data)
# print("pdb_no_amplifed_bp_GCN_512_gap-gmp_0.2_8.0.pkl = ",loaded_data1)


