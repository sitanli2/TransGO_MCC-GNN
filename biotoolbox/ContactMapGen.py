# import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
import os
# from contact_map_builder import ContactMapContainer,DistanceMapBuilder

# mapbuilder = DistanceMapBuilder
# strtest = "MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNTPYCFDQLRRRFGDVFSLQLAWTPVVVLNGLAAVREAMVTRGEDTADRPPAPIYQVLGFGPRSQGVILSRYGPAWREQRRFSVSTLRNLGLGKKSLEQWVTEEAACLCAAFADQAGRPFRPNGLLDKAVSNVIASLTCGRRFEYDDPRFLRLLDLAQEGLKEESGFLREVLNAVPVLPHIPALAGKVLRFQKAFLTQLDELLTEHRMTWDPAQPPRDLTEAFLAKKEKAKGSPESSFNDENLRIVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQRGRRVSPGCPIVGTHVCPVRVQQEIDDVIGQVRRPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSRDIEVQGFRIPKGTTLITNLSSVLKDEAVWKKPFRFHPEHFLDAQGHFVKPEAFLPFSAGRRACLGEPLARMELFLFFTSLLQHFSFSVAAGQPRPSHSRVVSFLVTPSPYELCAVPR"
# print(len(strtest)) #验证一下是不是515个氨基酸残基

def ContactMap_densities_Switching(changeType='densification',K_de=1,K_sp=3,initial_filename = '1A0A-A_P07270.npy'): #邻接矩阵稀疏化/密集化
    if (changeType == 'densification'): #邻接矩阵密集化
        """
        增加连接密度（接触图密集化）：若与残基A有接触关系的氨基酸残基集和为{A},与残基B有接触关系的氨基酸残基集和为{B},
        若A交B大于K，且A与B没有接触关系，此时视A与B相互接触。
        具体来说，就是两个残基有足够多的共同邻居，但是他们本身并不接触，此时也可以视他们为接触。
        """
        numpyArray = "G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/contact_map_dense_pdb_8.0/" + initial_filename
        data = np.load(numpyArray)

        # testlist = [[0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3],
        #              [0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6]]
        # data = np.array(testlist)

        data_list = data.tolist() #将numpy数组转换为列表（source,target形式）
        source_list = data_list[0]
        target_list = data_list[1]
        Amino_num = max(source_list)
        i = 0
        j = 0
        print(source_list,"\n",target_list)
        mapping_dict = {}

        while (i < len(source_list)):
            list1 = []  # 用来存储此时的key下的所有target
            key = source_list[i]  # 给一个初始值i.e列表的第一个元素
            while(j <= len(source_list)): #这个=号很关键，最后一个key不能漏了
                if (j == len(source_list)):  # 此时指针j已经指向最后一个元素了，再加1就要跳出循环了
                    mapping_dict[key] = list1
                    i = j
                    break
                elif(source_list[j] == key): #说明此时指针只想的元素任然属于key
                    list1.append(target_list[j])
                    j += 1
                else: #此时指针j指向了一个新的target
                    mapping_dict[key] = list1
                    i = j #此时指针i直接跳到j的位置
                    break
        print(mapping_dict)

        '''计算关系集和的平均长度'''
        # count = 0
        # total = 0
        # for key,value in mapping_dict.items():
        #     total += len(value)-1 #记得减去和自己的关系，算有效的总关系个数total
        #     count += 1
        # mean = total / count
        # print(count,total,mean)


        '''计算mappingdict 中每一个氨基酸和其它氨基酸的置信度（相似度）'''
        similarity_list = []
        for key1 in mapping_dict.keys():
            mother_list = mapping_dict[key1]
            #score_list = []
            if (len(mapping_dict[key1]) == 1):  # 自闭症，只和自己有关系
                # score_list.append(key1)
                similarity_list.append(key1)
                continue
            for key2 in mapping_dict.keys():
                if(len(mapping_dict[key2]) == 1):# 自闭症，只和自己有关系
                    continue
                if(key1 == key2):# 和自己的相似度不用比
                    continue
                if(key2 in mapping_dict[key1]):# 就像比血缘关系一样，已经存在（母子）关系了，不用比了
                    continue
                son_list = mapping_dict[key2] #key1和key2此时无直接关系，此时来比他们的相似度（看他们有的共同亲戚，在所有亲戚中所占的比例）
                common_elements = set(mother_list) & set(son_list)
                count = len(common_elements) # 找到他们之间的公共亲戚个数
                # w_mother = len(mother_list)/mean #计算权值，即若与其拥有关系的节点很少，score也会做适当修正
                # w_son = len(son_list)/mean
                # score_mother = w_mother * count/(len(mother_list)-1) #和mother_list有公共亲戚的个数除以mother_list亲戚的个数
                # score_son = w_son * count/(len(son_list)-1)
                score_mother = count / (len(mother_list) - 1)  # 和mother_list有公共亲戚的个数除以mother_list亲戚的个数
                score_son = count / (len(son_list) - 1)
                score = score_mother+score_son
                # score_list.append([key1,key2,score])
                similarity_list.append([key1,key2,score])
        # print(similarity_list) #输出包含计算得分的残基对：[source,target,score]

        # 根据第三个元素(相关性得分)进行降序排序
        sorted_list = sorted(similarity_list, key=lambda x: x[2], reverse=True)

        # 输出排序后的列表
        print(sorted_list)
        choose_num = 0
        final_list = []
        for hide_edge in sorted_list:
            if (choose_num >= K_de):#选择需要找到几条最相关的边,一般来说选5*2条比较好,*2是因为边有从source->target也有从target->source
             break
            edge_index = hide_edge[:2] #前两个元素就是两个有隐藏连接关系的节点,如:[34, 37, 1.4666666666666668] ,1.4666是他们的相关性系数
            if (edge_index not in final_list and edge_index[::-1] not in final_list): #就是看[node1,node2]和[node2,node1]其中有一个存过了不需要再存了
                final_list.append(edge_index)
                choose_num += 1
        print(final_list)
        S = source_list #初始source list
        T = target_list
        for edge_item in final_list: #循环递归
           print("hidden edge:",edge_item)
           S,T = insert_hidden_edge(source_list = S, target_list = T, hidden_edge = edge_item)
           S,T = insert_hidden_edge(source_list = S, target_list = T, hidden_edge = edge_item[::-1]) #source，target互换
        print(S,'\n',T)
        adjacency_matrix = [S, T]
        ContactMap = np.array(adjacency_matrix)
        print("符合GCN输入的adjacency_matrix：\n", ContactMap)

        # save_name = 'testfile'
        # contactmap_savepath = "G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/contact_map_dense_densification_8.0/"
        # filepath = contactmap_savepath + save_name
        # np.save(filepath, ContactMap)
        return ContactMap

    elif (changeType == 'sparsification'):
        '''减小连接密度（接触图稀疏化）：K近邻聚类，只保留接触距离最近的K个残基接触关系。'''
        print("K近邻聚类，只保留接触距离最近的K个残基接触关系")



def insert_hidden_edge(source_list=[0, 0, 0, 1, 1, 1, 1, 1, 1],target_list=[0, 1, 2, 0, 1, 2, 3, 4, 6],hidden_edge=[1,5]):
    '''用来按顺序insert元素进source,target
    有一个int元素组成的二维列表：list=[[source_list],[target_list]],如：
    [[0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
     3, 4, 4, 4, 4, 4, 4, 4,5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7]，
     [0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5,
     6, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 7, 8, 9, 10, 4, 5, 6, 7]]
     其中两个内嵌列表list1,list2中的元素都是一一对应的，source_list是有序列表，target_list是局部有序列表
    现在要将新的边索引元素[source1,target1]插入列表list中的两个内嵌列表中而不改变有序性，具体来说，source1找到有序列表list1应该插入的位置，随后再找target1的局部有序性，
    最后插入source1,target1进source_list,target_list
    '''
    # source_list = [0, 0, 0, 1, 1, 1, 1, 1, 1] #用来测试
    # target_list = [0, 1, 2, 0, 1, 2, 3, 4, 5] #用来测试
    # hidden_edge = [0,3] #用来测试

    source = hidden_edge[0]
    target = hidden_edge[1]
    index = 0
    if_append = 0

    while(index <= len(source_list)-1):
        if(source_list[index] == source): #找到有序列表source_list中source的首位（后面再从这个位置找target_list中的局部有序性）
            break
        index += 1

    while(source_list[index] == source): #找局部有序性
        if(target_list[index]>target): #找到target要插入的位置
            break
        if(index == len(source_list)-1): #可能index是在列表的队尾
            if_append = 1 #直接append了
            break
        # if(source_list[index+1] != source): #已经到了局部有序性的边缘，再往后走就是下一个source了
        index += 1

    if(if_append == 1):
        source_list.append(source)
        target_list.append(target)
    else:
        #此时要在index位置插入source,target
        source_list.insert(index,source)
        target_list.insert(index,target)

    '''
    target_list = [1, 2, 5, 8, 9]
    # 在索引位置 3 插入新元素 7
    target_list.insert(3, 7) #my_list.insert(index, element)：将 element 插入到index=3位置，即target_list[3]=8 ，原位置及之后的元素(8,9)会向后移动。
    target_list = [1, 2, 5, 7, 8, 9]
    '''
    # print(source_list,'\n',target_list)
    return source_list,target_list





def contact_mapping_AlphaFold(folder_path = 'G:/fkfk/AlphaFold_refine/biotoolbox/AlphaFoldfile_Dataset_Ultimate',threshold=8.,GCN_transformer=1):
    # folder_path为存放从AlphaFold下载的蛋白质结构的文件夹路径
    count = 1
    error_PDB_file_list=[] # 有些时候，PDB文件会有损坏，比如两个坐标间空格丢失：15.724-102.443 ,这个时候，提取坐标切片会出错，我们需要跳过这些坏文件
    for filename in os.listdir(folder_path):
        new_filename = filename[:-4] #去掉后缀 .pdb
        list1 = new_filename.split("_")
        print(list1,new_filename,f"filecount={count}")
        PDBID = str(list1[0])
        UPID = str(list1[1])
        start_point = int(list1[2])
        end_point = int(list1[3])
        save_name = PDBID+"_"+ UPID + '.npy'#即最终生成的contactmap 的文件名 i.e: 1A02-J_P05412.npy
        f = open("G:/fkfk/AlphaFold_refine/biotoolbox/AlphaFoldfile_Dataset_Ultimate/" + str(filename), "r") #读取该PDB文件
        D = [] #用来存放所有Ca原子的三维坐标
        """提取Ca原子坐标"""
        try:
            for a in f.readlines():  # 遍历该PDB文件的每一行
                b = a.split()  # 对这一行切片
                if (b[0] == "ATOM" or b[0] == "HETATM"):  # 找到 "ATOM" 打头的行，和"HETATM"的行，这里一定不要漏了HETATM ！！！
                    """
                    "ATOM" 打头的行用于表示：蛋白质或核酸的标准残基中的原子。
                    "HETATM"打头的行用于表示：用于表示：非标准残基或配体中的原子。
                    """
                    if(int(b[5])<start_point): #从蛋白质链的起始氨基酸位置开始提取，b[5] 指对应的行的原子属于蛋白质链中的第几个氨基酸（position）
                        continue

                    if(int(b[5])>end_point): #氨基酸残基终止位置
                        break

                    if (b[0] == "TER"):  # "TER" 打头的行表示一个链的结束。它的主要作用是指示某个链的终止位置，特别是在同一结构中有多个链时。
                        break


                    if b[2] == "CA":  # 找到 "CA" 原子的那一行
                        """
                        # 通过循环提取出每一个 CA 原子的空间三维坐标,但是有的时候，PDB文件中易出现两列并在一起的状况，此时切片序号会少一项数，导致定位到后一项'1.00170.93'
                        # 此时6，7，8需要换成5，6，7，所以坐标的起始值不应该固定为5，要我们自己求出来，因为三个坐标都是连在一起的，所以知道第一个的位置，就知道其他两个
                        position = 0 #用来定位第一个坐标的位置
                        while(position <= len(b)):
                            if b[int(position)].count('.') == 1 and b[int(position)].replace('.', '').isdigit(): #从b[]的第一个元素 'ATOM'开始往下找第一个满足这个要求的就是坐标起始点
                            # 只有一个小数点，并且其他字符都是数字（这样才是坐标），才能转浮点数
                                break
                            else:
                                position += 1
                        D.append((float(b[int(position)]),float(b[int(position)+1]),float(b[int(position)+2])))
                        """

                        D.append((float(b[6]), float(b[7]), float(b[8])))
                        # amino_list.append(b[5])  # b[5] 指对应的行的原子属于蛋白质链中的第几个氨基酸（position）
                        # amino += 1
        except:
            print(f"文件损坏，提取{filename}的Ca原子坐标时有误")
            error_PDB_file_list.append(str(filename))
            continue #跳过该损坏文件，解析下一个pdb文件
        """
           接下来通过两重循环计算每个氨基酸与任意一个氨基酸CA原子的距离
        """
        distance1 = []
        for b in range(len(D)):
            distance2 = []
            for c in range(len(D)):
                distance = ((D[b][0] - D[c][0]) ** 2 + (D[b][1] - D[c][1]) ** 2 + (
                            D[b][2] - D[c][2]) ** 2) ** 0.5  # 其实就是欧氏距离计算公式
                if distance == 0.: #自己和自己肯定接触
                    distance2.append(0)
                elif (distance <= threshold and distance != 0):
                    distance2.append(
                        round(float(distance), 3))  # round(float(distance), 3) 这段代码是将一个浮点数 distance 四舍五入保留三位小数
                elif (distance > threshold):  # 检查distance是否大于threshold，大于接触阈值则不接触，至为零
                    distance2.append(0)
            distance1.append(distance2)
        # print(distance1)
        if (GCN_transformer == 0):
            """
            将矩阵(n*n的列表)直接转换为numpy形式（全连接图）
            """
            ContactMap = np.array(distance1)
            print("不符合GCN输入的完整adjacency_matrix,i'e:全连接图：", ContactMap)
            filepath = "./Test_ContactMapFiles/" + PDBID + "_" + str(threshold) + "A-intact.npy"
            np.save(filepath, ContactMap)
            # 加载保存的 .npy 文件，以验证保存是否成功
            loaded_data = np.load(filepath)
            print("Loaded data:\n", loaded_data)
        else:
            """
            将矩阵转换为source,target 形式，作为GCN的邻接矩阵输入
            """
            distance1 = np.array(distance1)  # 转为numpy矩阵
            print("value of the generated contact map：\n", distance1, "\n shape of the map：\n", np.shape(distance1))
            Source = []
            Target = []

            """遍历接触图，找到存在接触的节点对"""
            for i in range(distance1.shape[0]):  # distance1.shape[0]就是取shape的第一个维度，即一共拥有氨基酸的个数
                for j in range(distance1.shape[1]):  # 通过两个for循环遍历整个 n*n的矩阵
                    if (distance1[i, j] != 0. or i == j):  # 判断第i行第j列是否为0，不为零则说明第i个氨基酸和第j个氨基酸残基互相接触,还有自己和自己也是接触的
                        # if(i<79 and j<79): #要小于链长，应为只从PDB文件中取其中的A/B链的一种，即结构的一部分。
                        Source.append(i)
                        Target.append(j)

            adjacency_matrix = [Source, Target]
            ContactMap = np.array(adjacency_matrix)
            print("符合GCN输入的adjacency_matrix：", ContactMap)
            contactmap_savepath = "G:/fkfk/AlphaFold_refine/UltraDataset/testData3_290entries/contact_map_dense_pdb_8.0/"
            filepath = contactmap_savepath + save_name
            np.save(filepath, ContactMap)
            count += 1
            # 加载保存的 .npy 文件，以验证保存是否成功
            # loaded_data = np.load(filepath)
            # print("Loaded data:\n", loaded_data)
    print("损坏文件列表(缺失空格，无法解析提取Ca原子坐标)：",error_PDB_file_list)
        # break #纯测试，第一轮结束后退出





def contact_mapping_PDB(PDBID='3UDW-C_P15151',threshold=8.,GCN_transformer=1,savepath = "G:/fkfk/AlphaFold_refine/biotoolbox/Contactmapfile_Dataset_Ultimate/"): #设定GCN_Transformer为0 ，即不做GCN数据转换（N*N ---> 2*r），生成完整的N*N接触图（方便我们可视化接触图）
    #测试 UPID = "A0A087X1C5" ,threshold为接触阈值默认为8, GCN_transformer代表是否要将接触图转换为source->target的形式
    #测试 3AQF-A_O60895 ,即UPID=O60895. PDBID=3AQF-A
    # f_test = open("./TestPDBDataset/"+str(PDBID)+".pdb","r")#PDB 文件
    f = open("G:/fkfk/AlphaFold_refine/biotoolbox/PDBfile_Dataset_Ultimate/"+str(PDBID)+".pdb","r")
    # f = open("./TestPDBDataset/"+"pdb"+str(PDBID)+".ent","r") #ent 文件

    D = []
    nextchain = False # 如果要下一条蛋白质链，则将此parameter设置为True
    amino = 56  # 该条链的起始氨基酸在全链中的位置
    amino_list = []  # 用来存放已经提取出Ca原子坐标的氨基酸残基位置
    for a in f.readlines(): #遍历PDB文件的每一行
        b = a.split() #对这一行切片
        if (b[0] == "TER"): #"TER" 打头的行表示一个链的结束。它的主要作用是指示某个链的终止位置，特别是在同一结构中有多个链时。
            break

        if (b[0] == "ATOM" or b[0] == "HETATM"): # 找到 "ATOM" 打头的行，和"HETATM"的行，这里一定不要漏了HETATM ！！！
            """
            "ATOM" 打头的行用于表示：蛋白质或核酸的标准残基中的原子。
            "HETATM"打头的行用于表示：用于表示：非标准残基或配体中的原子。
            """
            if b[2] == "CA": # 找到 "CA" 原子的那一行
                """
                # 通过循环提取出每一个 CA 原子的空间三维坐标,但是有的时候，PDB文件中易出现两列并在一起的状况，此时切片序号会少一项数，导致定位到后一项'1.00170.93'
                # 此时6，7，8需要换成5，6，7，所以坐标的起始值不应该固定为5，要我们自己求出来，因为三个坐标都是连在一起的，所以知道第一个的位置，就知道其他两个
                position = 0 #用来定位第一个坐标的位置
                while(position <= len(b)):
                    if b[int(position)].count('.') == 1 and b[int(position)].replace('.', '').isdigit(): #从b[]的第一个元素 'ATOM'开始往下找第一个满足这个要求的就是坐标起始点
                    # 只有一个小数点，并且其他字符都是数字（这样才是坐标），才能转浮点数
                        break
                    else:
                        position += 1
                D.append((float(b[int(position)]),float(b[int(position)+1]),float(b[int(position)+2])))
                """

                D.append((float(b[6]), float(b[7]), float(b[8])))
                amino_list.append(b[5]) #b[5] 指对应的行的原子属于蛋白质链中的第几个氨基酸（position）
                amino += 1

    if (nextchain):
        D_nextchain = []
        flag = 0
        for a in f.readlines(): #遍历PDB文件的每一行
            b = a.split() #对这一行切片
            if (b[0] == "TER"): #"TER" 打头的行表示一个链的结束。它的主要作用是指示某个链的终止位置，特别是在同一结构中有多个链时。
                flag = 1
                continue
            if (b[0] == "TER"):
                break
            if (b[0] == "ATOM" or b[0] == "HETATM"): # 找到 "ATOM" 打头的行
                if b[2] == "CA": # 找到 "CA" 原子的那一行
                    if(flag == 0):
                        D.append((float(b[6]),float(b[7]),float(b[8]))) # 通过循环提取出每一个 CA 原子的空间三维坐标
                        #print(len(D))
                    if(flag == 1):
                        D_nextchain.append((float(b[6]),float(b[7]),float(b[8])))
        chain1_len = len(D)
        chain2_len = len(D_nextchain)
    """
    接下来通过两重循环计算每个氨基酸与任意一个氨基酸CA原子的距离
    """

    distance1 = []
    for b in range(len(D)):
        distance2 = []
        for c in range(len(D)):
            distance = ((D[b][0] - D[c][0]) ** 2 + (D[b][1] - D[c][1]) ** 2 + (D[b][2] - D[c][2]) ** 2) ** 0.5 #其实就是欧氏距离计算公式
            if distance == 0.:
                distance2.append(0)
            elif (distance <= threshold and distance != 0):
                distance2.append(round(float(distance),3)) #round(float(distance), 3) 这段代码是将一个浮点数 distance 四舍五入保留三位小数
            elif (distance > threshold): # 检查distance是否大于threshold
                distance2.append(0)
        distance1.append(distance2)
    #print(distance1)
    if(GCN_transformer == 0):
        """
        将矩阵转换为numpy形式
        """
        ContactMap = np.array(distance1)
        print("不符合GCN输入的完整adjacency_matrix：", ContactMap)
        filepath = "./Test_ContactMapFiles/" + PDBID +"_"+ str(threshold) +"A-intact.npy"
        np.save(filepath, ContactMap)
        # 加载保存的 .npy 文件，以验证保存是否成功
        loaded_data = np.load(filepath)
        print("Loaded data:\n", loaded_data)
    else:
        # print("value of the generated contact map：\n",distance1,"\n shape of the map：\n",np.shape(distance1))
        # 经打印后可知，这个contact map 是一个515*515 即n*n的矩阵，我还需要将他转换为符合GCN输入的[Source,Target] 即2*n的形式作为整个GCN的edge_index
        """
        将矩阵转换为numpy形式
        """
        distance1 = np.array(distance1)  # 转为numpy矩阵
        print("value of the generated contact map：\n", distance1, "\n shape of the map：\n", np.shape(distance1))
        Source = []
        Target = []

        """遍历接触图，找到存在接触的节点对"""
        for i in range(distance1.shape[0]): #distance1.shape[0]就是取shape的第一个维度，即一共拥有氨基酸的个数
            for j in range(distance1.shape[1]): #通过两个for循环遍历整个 n*n的矩阵
                if (distance1[i,j] != 0. or i == j): # 判断第i行第j列是否为0，不为零则说明第i个氨基酸和第j个氨基酸残基互相接触,还有自己和自己也是接触的
                    # if(i<79 and j<79): #要小于链长，应为只从PDB文件中取其中的A/B链的一种，即结构的一部分。
                        Source.append(i)
                        Target.append(j)

        adjacency_matrix = [Source,Target]
        ContactMap = np.array(adjacency_matrix)
        print("符合GCN输入的adjacency_matrix：",ContactMap)
        saved_filename = PDBID+ ".npy" #这里其实就是文件保存的名字："3UDW-C_P15151.npy"
        filepath = savepath + saved_filename
        np.save(filepath,ContactMap)

        # 加载保存的 .npy 文件，以验证保存是否成功
        loaded_data = np.load(filepath)
        print("Loaded data:\n",loaded_data)



def UseBiopython():
    """使用PDBparser读取和解析PDB文件"""
    pdb_code = "A0A087X1C5"
    parser = PDBParser()
    structure = parser.get_structure(pdb_code, "./TestPDBDataset/A0A087X1C5.pdb")

    # 计算残基之间的距离
    def calc_residue_dist(residue_one, residue_two):
        diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

    def calc_dist_matrix(chain_one, chain_two):
        dist_matrix = np.zeros((len(chain_one), len(chain_two)), dtype=np.float)
        for i, residue_one in enumerate(chain_one):
            for j, residue_two in enumerate(chain_two):
                dist_matrix[i, j] = calc_residue_dist(residue_one, residue_two)
        return dist_matrix

    # 构建接触图
    chain_one = structure[0]["A"]
    chain_two = structure[0]["B"]
    dist_matrix = calc_dist_matrix(chain_one, chain_two)
    contact_map = dist_matrix < 12.0  # 使用阈值确定接触

    # 打印最小和最大距离
    print("Minimum distance:", np.min(dist_matrix))
    print("Maximum distance:", np.max(dist_matrix))
    print("contact map:",contact_map)


    """通过networkx 和 matplotlib进行可视化"""
    # G = nx.Graph()
    #
    # nx.from_numpy_matrix(distance1)
    # nx.draw(G)
    # plt.show()
    # nx.betweenness_centrality(G)
    # nx.closeness_centrality(G)
    # nx.degree_centrality(G)
    # nx.clustering(G)

def verify():
    array1 = np.load("G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/contact_map_dense_pdb_8.0/3AQF-A_O60895.npy")
    array2 = np.load("G:/fkfk/Struct2Go-main/data_collect/contact_map_dense_pdb_8.0/3AQF-A_O60895.npy")
    # 比较两个数组是否一致
    are_equal = np.array_equal(array1, array2)
    if are_equal:
        print("两个 .npy 文件的内容一致。")
    else:
        print("两个 .npy 文件的内容不一致。")

def list_files_in_directory(directory):
    # 遍历指定目录下的所有文件
    file_name_list = []
    for filename in os.listdir(directory):
        # 获取文件的完整路径
        file_path = os.path.join(directory, filename)
        # 检查是否是文件
        if os.path.isfile(file_path):
            file_name_list.append(filename)
            print(filename)
    print('total file number=',len(file_name_list))
    return file_name_list


if __name__ == '__main__':

    """Biopython 工具包测试"""
    # UseBiopython()

    """自己写的接触图生成函数（PDB,AlphaFold）"""
    # contact_mapping(UPID="3aqf",threshold=8.,GCN_transformer=1)# 设定阈值为6，8，10，12.  设定GCN_Transformer为0 ，即不做GCN数据转换（N*N ---> 2*r），生成完整的N*N接触图（方便我们可视化接触图）
    # contact_mapping_AlphaFold(folder_path = 'G:/fkfk/AlphaFold_refine/biotoolbox/AlphaFoldfile_Dataset_Ultimate',threshold=8.,GCN_transformer=1)


    """插入隐藏边进入[[source],[target]]"""
    # insert_hidden_edge()

    """生成添加隐藏边（或删去冗余边）后的接触图"""
    # file_name_list = list_files_in_directory("G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/contact_map_dense_pdb_8.0")
    # error_file_list = []
    # for file_name in file_name_list:
    #     """接触图隐藏边密集化扩充 or 稀疏化"""
    #     try:
    #         ContactMap = ContactMap_densities_Switching(initial_filename=str(file_name),K_de=1)
    #
    #         """生成densification后的contact map 并存入新的文件夹"""
    #         save_name = file_name
    #         contactmap_savepath = "G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/contact_map_dense_densification_8.0/"
    #         filepath = contactmap_savepath + save_name
    #         np.save(filepath, ContactMap)
    #     except: #小部分接触图再寻找隐藏边时会有索引问题（比如有氨基酸残基只和自己接触的情况），此时出错的接触图咱们先保存他的file_name
    #         error_file_list.append(file_name)
    #         print('error file=',error_file_list)
    #
    # print('error file=', error_file_list)
    # file_name_list = list_files_in_directory("G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/contact_map_dense_densification_8.0")

    """验证两个[[source],[target]]形式的接触图是否相同"""
    verify()