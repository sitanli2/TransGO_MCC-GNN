from pymongo import MongoClient
import pandas as pd
import numpy as np
import os
import requests
from AlphaFold_refine.biotoolbox.ContactMapGen import contact_mapping_PDB
# import tensorflow as tf


def MongoConnection(collection_name = 'ProGO_Target_Info'):
    mongoClient = MongoClient('mongodb://localhost:27017/')
    database = mongoClient.get_database('dtmr_dev')

    # collection = dtmr_dev['uniprot_Target_Info_test']

    collection = database.get_collection(collection_name)
    return collection

    # # 建立mongodb连接
    # mongoClient = MongoClient('mongodb://localhost:27017/')
    # # mongoClient = MongoClient('mongodb://172.10.10.9:27017/')
    # dtmr_dev = mongoClient.get_database('dtmr_dev')
    # # dtmr_dev.authenticate('dtmr_dev', 'dtmr_dev') #密码验证
    #
    # dev = dtmr_dev.get_collection('ProGO_Target_Info') #uniprot_Target_Info_test
    # dev2 = dtmr_dev.get_collection('ProGO_skip_error_id_list')
    # devlog = dtmr_dev.get_collection('ProGO_Log')


def return_id_name_list(data):
    i = 0#行指针
    j = 0#列指针(这个csv读取后生成的dataframe文件只有两列，一列id一列name)
    uniprot_id_list = []
    uniprot_name_list = []
    hang_len = data.shape[0]
    lie_len = data.shape[1]
    while(i <= hang_len-1 and j <= lie_len-1):
        uniprot_id_list.append(data.loc[i][j])
        i += 1#换行
        j += 1#换到第二列得到name
        uniprot_name_list.append(data.loc)
        j = 0 #列换回到第一列id列
    return uniprot_id_list,uniprot_name_list

def Get_UPIDlist(file_path = 'G:/fkfk/AlphaFold_refine/UltraDataset/Database1+/IDlist.txt'):
    UPID_list = []
    with open(file_path, 'r') as file:
        # 读取第一行标题并忽略
        title = file.readline().strip()

        # 逐行读取剩下的内容
        for line in file:
            # 去除每行的前后空格并按空格分割成列表
            elements = line.strip()
            UPID_list.append(elements)
        # print(len(UPID_list),'\n',
        #       UPID_list)
    return UPID_list

def Gen_GO_labelTXT(UPID_list = []):
    """用来给MonGoDB中的每一个对应条目都加上mf,bp,cc"""
    # UPID_list = Get_UPIDlist(file_path ='/AlphaFold_refine/Protein_spider/UniProt_Scrapy/uniprot_id_expanded')
    UPID_list = TXT_to_List('G:/fkfk/AlphaFold_refine/Protein_spider/UniProt_Scrapy/DeepFRI_mapping_UPID')
    collection = MongoConnection(collection_name = 'ProGO_Target_Info')
    # data = pd.read_csv('G:/fkfk/AlphaFold_refine/Protein_spider/UniProt_Scrapy/uniprot_id.csv', engine='python')  # sep='/n',
    # uniprot_id_list, uniprot_name_list = return_id_name_list(data) #这个uniprot_id_list目前还不是完整的，要再加上No repeat UPID from Database1+
    testcount = 0

    FullSize = int(len(UPID_list))
    print("Mongo中现有的条目数=",FullSize)
    TestSize = int(input("输入你想要在Mongo中更新（加上mf,bp,cc）的条目数"))
    for UPID in UPID_list:
        errorcode = 0
        try:
            query_result = collection.find_one({"UPID": UPID})["Go_info"]  # 从MonGoDB中索引到一个条目的 Go_info,
            # 这里需要注意一个问题，有可能出现errorID的情况（网络接口返回文件损坏如UPID=A6H8Y1），详情见爬虫部分，
            # 我把这5个errorID单独存起来了，所以当遍历到它们的ID时，靠上面这行代码检索将会出错，应为他们在Mongo中为null
            # 所以在这里可以给这个查询加一个判错，当遍历到errorID（查询失败时），自动跳过。
            print("查询结果:", query_result, '\n', type(query_result))
        except Exception as error:
            print('错误类型是', error.__class__.__name__)
            print('错误明细是', error)
            errorcode = 1


        MFlist = []
        BPlist = []
        CClist = []
        if (len(query_result) > 0 and errorcode == 0):  # 不为空列表，即拥有Go_info（一般来说肯定都是有的，除非爬虫出问题了导致上面的查询出错）
            for item in query_result:
                Golabel = item["id"]
                property = item["properties"]
                for prop in property:
                    if (prop["key"] == "GoTerm"):
                        if (str(prop["value"])[0] == "C"):  # 说明这个Golabel是一个CC 标签
                            CClist.append(Golabel)
                        elif (str(prop["value"])[0] == "F"):  # 说明这个Golabel是一个MF 标签
                            MFlist.append(Golabel)
                        elif (str(prop["value"])[0] == "P"):  # 说明这个Golabel是一个BP 标签
                            BPlist.append(Golabel)
                    break  # 这个循环主要是为了遍历整个 properties 列表，拿到了我们想要的GoTerm后就要跳出这个循环
        else:
            print("UPID=" + UPID + "的MonGoData没有Go标签条目")
            pass #存入空列表 mf,bp,cc

        # 通过上面的一系列操作后，我们拿到了这条UPID对应的 所有Go_label,并将它们分别存进了对应的三个category的列表中
        record = {}
        record["mf"] = MFlist
        record["bp"] = BPlist
        record["cc"] = CClist

        # 更新文档
        update_result = collection.update_one({"UPID": UPID}, {"$set": record})  # 给每个条目都加上record

        # 打印更新后的文档
        if update_result.modified_count > 0:
            print("更新成功,更新条目：", record)

        else:
            print("未找到对应的文档")

        """
        我们这样做的目的是为了通过下面的Generate_Go_label()函数生成最后的Go_label.txt 文件，element如下所示：

        >1E0T-A_P0AD61
        mf:GO:0016772	GO:0016773	GO:0016301	GO:0097367	GO:0042802	GO:0032559	GO:0000287	GO:0017076	GO:0032555	GO:0005524	GO:0035639	GO:0030554	GO:0032553	
        bp:GO:0006950	GO:0009135	GO:0046939	GO:0016052	GO:0019637	GO:0009185	GO:0051259	GO:0006163	GO:0065003	GO:0051262	GO:0019693	GO:0009179	GO:0043436	GO:0005975	GO:0006165	GO:0016310	GO:0006757	GO:0032787	GO:0009117	GO:0071840	GO:0006091	GO:0022607	GO:1901575	GO:0009056	GO:0006753	GO:0006082	GO:0009266	GO:0055086	GO:0009259	GO:0009628	GO:0046031	GO:0044085	GO:0016043	GO:0009408	GO:0006096	GO:0006796	GO:0051260	GO:0043933	GO:0009150	GO:0009132	GO:0006793	GO:0072521	GO:1901135	GO:0019752	GO:0051289	GO:0006090	GO:0046034	
        cc:GO:0005829
        以测试数据 即：UPID=A0A087X1C5 为例：
        {'mf': ['GO:0070330', 'GO:0020037', 'GO:0005506', 'GO:0016712'], 'bp': ['GO:0019369', 'GO:0042178', 'GO:0006805'], 'cc': ['GO:0005737', 'GO:0043231', 'GO:0016020', 'GO:0005739']}

        于是它在Go_label.txt中：
        >A0A087X1C5 （AF-A0A087X1C5-F1 这个ID在MongoDB中的Structure_AlphaFold条目下可以找到,可以通过这个ID爬取AlphaFold的PDB文件）
        mf:GO:0070330   GO:0020037  GO:0005506  GO:0016712
        bp:GO:0019369   GO:0042178  GO:0006805
        cc:GO:0005737   GO:0043231  GO:0016020  GO:0005739
        """

        testcount += 1  # 用来标记测试数据的条数
        if (testcount == TestSize):
            break

# def Gen_GO_labelTXT_usecursor():
#     TestSize = int(input("输入你想要在Mongo中更新（加上mf,bp,cc）的条目数"))
#     collection = MongoConnection()
#     cursor = collection.find()  # 获取一个游标（cursor），然后逐个迭代游标以按顺序遍历MongoDB中的每个条目访问每个文档。
#     for document in cursor:



def GetFullIDlist():
    """
    为了避免重复 根据本地的Database1+ 的 pdb_mapping_seq.txt文件获取所有的七千多条UPID 和 PID
    :return:FullUPIDlist,FullPIDlist
    """
    path = "/AlphaFold_refine/biLSTM+GCN/data_collect/pdb_mapping_seq.txt"
    FullUPIDlist = []
    FullPIDlist = []
    with open(path,"r") as file:
        for line in file:
            if '>' in line:
                pid_chain_uid = line[1:].strip('\n').split('\t')[0]
                PID = pid_chain_uid.split("_")[0]
                UPID = pid_chain_uid.split("_")[1]
                FullUPIDlist.append(UPID)
                FullPIDlist.append(PID)
                # print("test data UPID =",UPId,"test data PID =",PID)
                # break
        #print("test data UPID len=", len(FullUPIDlist), "test data PID len=", len(FullPIDlist))
    return FullUPIDlist,FullPIDlist


def getAlphaEntryIDFromMongo(UPID):
    collection = MongoConnection(collection_name = 'ProGO_Target_Info')
    try:
        query_result = collection.find_one({"UPID": UPID})["Structure_AlphaFold"]["IDENTIFIER_entryId"]  # 从MonGoDB中通过UPID索引到其对应条目的 Structure_AlphaFold 嵌套字典
        # print(UPID+"的IDENTIFIER_entryId 查询结果:", query_result,'\n返回类型：', type(query_result))
        entryID = query_result
    except:
        print("查询不到与UPID=",UPID,"对应的AlphaFold Entry条目，没有与之对应的AlphaFold预测结构")
        entryID = 'null'
    return entryID

def getPDBEntryIDFromMongo(UPID):
    collection = MongoConnection(collection_name = 'ProGO_Target_Info')
    try:
        query_result = collection.find_one({"UPID": UPID})["Structure_PDB"]["IDENTIFIER_entryId"]  # 从MonGoDB中通过UPID索引到其对应条目的 Structure_AlphaFold 嵌套字典
        # print(UPID+"的IDENTIFIER_entryId 查询结果:", query_result,'\n返回类型：', type(query_result))
        entryID = query_result
    except:
        print("查询不到与UPID=",UPID,"对应的PDB Entry条目，没有与之对应的PDB精确测定结构")
        entryID = 'null'
    return entryID


def LabelhaveMf_andAll_InProGO(filterNum=100,mf_func_list=[],bp_func_list=[],cc_func_list=[]):
    """
    此函数用来访问Mongo数据库中爬虫爬取的数据，将他们对应的6000条数据的UPID保存下来生成 Intact_label.txt
    这6000条数据应该满足以下筛选条件：
    - 通过 GetFullIDlist() 获得Database1 的数据UPID,Database2筛选的数据不应该与Database1的数据发生重复
    - 由于ESM-1b 能处理的氨基酸序列最大长度为1022个氨基酸+首token+尾token  所以所有蛋白质序列长度超过1022的都不行，否则就需要：（为了省事，我就直接排除序列长度超过1022的蛋白质了）
        ESM - 1b 模型在不超过1024个残基限制（这1024个residue还包括了起始和终止的token,所以真正的训练残基数不能超过1022个）的条件下进行训练
        所以对于长度大于1022的蛋白质序列，比如说序列长度为n，我们将序列分为 n / 1022 个长度为1022的块，除去最后一个长度为 n % 1022的块（去除余数），每一个块（1022长的氨基酸序列）都生成一个嵌入（embedding），
        最后将所有块的嵌入连接（concate 拼接，连接）起来作为这个序列长度超过1022的蛋白质的嵌入。
    - 筛选出来的蛋白质条目需要有高质量的蛋白质空间结构（PDB结构）
    - 最后按照mf,bp,cc都要有的条件来筛选出最终的6000条数据作为Database2

    :return:
    """
    IntactLabelID_list = []
    collection = MongoConnection()
    cursor = collection.find() #获取一个游标（cursor），然后逐个迭代游标以按顺序遍历MongoDB中的每个条目访问每个文档。
    FullUPIDlist,FullPIDlist = GetFullIDlist()

    # testcount = int(input("你需要从数据库中筛选出几条ID呢？"))
    accomplish = 0 #用来标志是否完成了设定的目标,i.e:testcount == len(IntactLabelID_list)
    for document in cursor:
        print(document,"\n",type(document))
        if (document["UPID"][0] != 'P'): #UPID稍微限制一下，首字母为"P"
            continue
        if (document["UPID"] in FullUPIDlist):
            print("该protein条目已经在Database1中存在了，跳过")
            continue
        errorID = ['A6H8Y1','P08579','Q06190','Q969Q6','Q99437','Q9Y5P8'] #这六个ID是在爬虫阶段产生错误的ID无法解析他们的数据接口返回值，需要后续单独用beautifulsoup或其他方法直接解析他们的网页HTML拿数据
        if (document["UPID"] in errorID):
            print("错误ID",document["UPID"])
            continue
        if ('inactive' in document) :
            # 数据库中还存在 未激活状态的数据比如：
            # {'_id': ObjectId('662deb4021712c0f197615b9'), 'UPID': 'Q08AF8', 'inactive': 'True', 'inactiveReasonType': ['DEMERGED'], 'bp': [], 'cc': [], 'mf': []}
            print("该数据在UniprotKB中未激活，inactiveReasonType = ",document['inactiveReasonType'])
            continue
        if (getAlphaEntryIDFromMongo(document["UPID"]) == 'null'): #在后续的测试中发现，有些蛋白质数据是没有AlphaFold预测结构的（比如UPID =A2VEC9、P22352..... ），这会导致数据缩水，所以要提前排除
            print("该UPID=",document["UPID"],"的蛋白质没有与之对应的AlphaFold结构")
            continue
            # break
        if (len(str(document["Sequence"])) >= 1022): #该蛋白质序列长度超过了Esm能处理的极限，跳过该蛋白
            continue

        if (len(document["Structure_PDB"]) == 0): #筛选出有高质量空间结构即PDBID的蛋白质条目
            continue
        PDBID = '' # PDB文件ID，如：6Y7F
        Chains_value = '' # 具体选的哪条蛋白质链，以及该蛋白质链在序列中所处的位置，如：A/B=1-281

        for PDB_structure in document["Structure_PDB"]: # 遍历所有PDB结构，选出我们要的
            try:
                if (PDB_structure["properties"][0]["value"] == "X-ray"): #筛选出有高质量空间结构即PDBID的蛋白质条目
                    PDBID = PDB_structure["id"]
                    Chains_value = PDB_structure["properties"][2]["value"]
                    break #符合这个判断条件就是找到了，记得break
            except:
                break

        if (PDBID == '' or Chains_value == ''): #没找到
            continue

        if ("mf" in document and "bp" in document and "cc" in document): #先检验一下，该条目下这三个文档是否成功通过上面的 Gen_GO_labelTXT()函数 添加
            #if ("Structure_AlphaFold" in document): #按理来说这些数据都会有AlphaFold entry的
            if(len(document["mf"]) >= 2 and len(document["bp"]) > 1 and len(document["cc"]) > 1):
                #这一步就是筛选看这个条目是否三个标签都存在，当然还可以更严格一点如要求三个标签的个数要大于某个阈值,我对mf标签比较严格只选取mf标签大于等于2个的条目
                common_mf = set(document["mf"]) & set(mf_func_list)
                common_bp = set(document["bp"]) & set(bp_func_list)
                common_cc = set(document["cc"]) & set(cc_func_list)
                if(bool(common_mf) and bool(common_bp) and bool(common_cc)):
                    # 将列表转换为集和，求集和是否有交集,即 看该UPID对应的蛋白质条目拥有的mf,bp,cc标签是否属于Database1+里面的Go_label,是的话，再加上这个UPID
                    print("该UPID对应的蛋白质条目与Database1+中的GoLabel 有mf交集：",common_mf,'\n',
                          "该UPID对应的蛋白质条目与Database1+中的GoLabel 有bp交集：",common_bp,'\n',
                          "该UPID对应的蛋白质条目与Database1+中的GoLabel 有cc交集：",common_cc,'\n',)
                    IntactID = str(PDBID) + "$" + str(Chains_value) + "$" + document["UPID"] # 生成完整的filteredID
                    IntactLabelID_list.append(IntactID)
        else:
            print("该条目还没通过Gen_GO_labelTXT()函数生成mf,bp,cc索引")
        if (filterNum == len(IntactLabelID_list)):
            accomplish = 1 #完成后置1
            break #这里我做了一个测试 没有break让他一直执行下去，最后得到的IntactLabelID_list 中有6459条数据，说明在爬下来的2万多条数据里面共有6459条数据满足上面的筛选条件
    if (accomplish == 1):
        print("value of the list:", IntactLabelID_list)
        print("成功筛选出指定",len(IntactLabelID_list),"条数据")
        return IntactLabelID_list
    else:
        print("value of the list:", IntactLabelID_list)
        print("爬虫基准数据库无法筛选出你想要的条数=",filterNum,"最多能筛选出的条数为=",len(IntactLabelID_list),"\n","我们需要继续扩充ProGO数据库！")
        return IntactLabelID_list


def List_to_TXT(testList,savePath): #换行写入
    with open(savePath,'w') as file:
        for item in testList:
            file.write(item + "\n")

def TXT_to_List(filePath): #换行读取成列表
    with open(filePath,'r') as file:
        IDlist = [line.strip() for line in file]
    # print(len(IDlist),IDlist)
    return IDlist




def GenIntact_labelTXT(IDlist,isAlphaFold=0):
    """
    生成以类似：1E0T-A_P0AD61 标签组成的文本文件
    isAlphaFold 为1 则在ProGO中查询P0AD61对应的AlphaFold entry ID
    为0 则询P0AD61对应的PDB entry ID
    :return:
    """
    savePath = "G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/pdb_mapping_seq_new.txt" #最终的生成路径
    # 将列表中的元素以（>4JOL-A_Q06455	EEMIDHRLTDREWAEEWKHLDHLLNCIMDMVEKTRRSLTVLRRCQEADREELNYWIRRYSDAE）的形式逐行写入到文本文件
    with open(savePath, 'w') as file:
        for item in IDlist:
            if(isAlphaFold == 1):
                AlphaFold_entryID = getAlphaEntryIDFromMongo(item)
                if(AlphaFold_entryID == 'null'):
                    continue #这条蛋白质没有合适的AlphaFold预测结构
                full_item = str(AlphaFold_entryID) +'_'+ str(item) #拼接起来，做为最终的ID
                file.write(full_item + '\n')
            else:
                PDB_entryID = ''



# def LabelhaveMf_andAll_Indisk(getnum):
#     """
#     此函数用来从data_collect 那个Go_label原件中筛选出所有有MF标签的IDlist 即 mfIDlist
#     随后再筛选出 BP,MF，CC 三个标签都有的IDlist 即 IntactLabelIDlist
#     :parameter:getnum IntactLabelIDlist需要几条 ，即拥有完整bp,cc,mf 的ID
#     :return: mfIDlist IntactLabelIDlist
#     """
#     IntactLabelIDlist = [] #用来存放mf,bp,cc三个标签都有的ID
#     mfIDlist = [] #用来存放有mf标签的这些ID
#     with open('../../data_collect/Go_label.txt', 'r') as f: #data_collect 下的Go_label 原件
#         for line in f:
#             if '>' in line: #说明这一行是ID 如 >4WKG-A_P77398，也作为新的一个开始需要把所有Flag置为0
#                 mfflag = 0
#                 bpflag = 0
#                 ccflag = 0
#                 pdb_chain_uid = line[1:].strip('\n')
#             elif 'mf' in line: #当前行是mf标签行
#                 if 'mf:' == line.strip('\n'): #说明这一整行 只有一个 ‘mf:’ ,即mf标签为空
#                     # NoMfIDlist.append(pdb_chain_uid)
#                     pass #当然也可以向上面一样存一下这个ID
#                 else:
#                     mfflag= 1 #将flag置为1
#                     mfIDlist.append(pdb_chain_uid)
#             elif 'bp' in line: #当前行是bp标签
#                 if 'bp:' != line.strip('\n'): #看这一行bp标签是否为空，不为空将对应的flag置为一
#                     bpflag = 1
#             elif 'cc' in line: #当前行是cc标签,经过几次循环后下一行有 '>' 这个时候就要判断一下几个flag了
#                 if 'cc:' != line.strip('\n'): #看这一行cc标签是否为空，不为空将对应的flag置为一
#                     ccflag = 1
#                 if (mfflag == 1 and bpflag == 1 and ccflag == 1): # 三个标签都有
#                     IntactLabelIDlist.append(pdb_chain_uid)
#             if(len(IntactLabelIDlist) >= getnum):
#                 break
#         return mfIDlist,IntactLabelIDlist
#
#
#
#
#
# def GenIntact_lebelTXT(IntactLabelIDlist):
#     filepath = "Intact_label.txt"
#     # 将列表中的元素逐行写入到文本文件
#     with open(filepath, 'w') as file:
#         for item in IntactLabelIDlist:
#             file.write(item + '\n')

def get_AllGOlabel(func = 'mf'):
    """
    这个函数用来得到 目标路径下的Go_label文件包含 Mf, BP, CC的具体信息,并以onehot的形式制作一个GO标签对照字典
    :param func: MF,BP,CC
    :return:
    """
    if func == 'mf':
        mf_dict = {}
        mf_func = []
        mf_label_dict = {}
        with open('/AlphaFold_refine/UltraDataset/Database1+/Go_label.txt', 'r') as f:
            for line in f:
                if '>' in line:
                    pdb_chain_uid = line[1:].strip('\n')
                elif 'mf' in line:
                    if 'mf:' == line.strip('\n'):
                        continue
                    else:
                        mf_func_list = line[3:].strip().split('\t')
                        mf_dict[pdb_chain_uid] = mf_func_list
        with open('/AlphaFold_refine/UltraDataset/Database1+/mf/mf_label.txt', 'r') as f:
            for line in f:
                line = line.strip('\n')
                mf_func.append(line)
        for i in mf_dict.keys():
            label = np.zeros(len(mf_func))
            for j in mf_dict[i]:
                if j in mf_func:
                    index = mf_func.index(j)
                    label[index] = 1
            mf_label_dict[i] = label
        return mf_label_dict,len(mf_func)
    elif func == 'bp':
        bp_dict = {}
        bp_func = []
        bp_label_dict = {}
        with open('/AlphaFold_refine/UltraDataset/Database1+/Go_label.txt', 'r') as f:
            for line in f:
                if '>' in line:
                    pdb_chain_uid = line[1:].strip('\n')
                elif 'bp' in line:
                    if 'bp:' == line.strip('\n'):
                        continue
                    else:
                        bp_func_list = line[3:].strip().split('\t')
                        bp_dict[pdb_chain_uid] = bp_func_list
        with open('/AlphaFold_refine/UltraDataset/Database1+/bp/bp_label.txt', 'r') as f:
            for line in f:
                line = line.strip('\n')
                bp_func.append(line)
        for i in bp_dict.keys():
            label = np.zeros(len(bp_func))
            for j in bp_dict[i]:
                if j in bp_func:
                    index = bp_func.index(j)
                    label[index] = 1
            bp_label_dict[i] = label
        return bp_label_dict,len(bp_func)
    elif func == 'cc':
        cc_dict = {}
        cc_func = []
        cc_label_dict = {}
        with open('/AlphaFold_refine/UltraDataset/Database1+/Go_label.txt', 'r') as f:
            for line in f:
                if '>' in line:
                    pdb_chain_uid = line[1:].strip('\n')
                elif 'cc' in line:
                    if 'cc:' == line.strip('\n'):
                        continue
                    else:
                        cc_func_list = line[3:].strip().split('\t')
                        cc_dict[pdb_chain_uid] = cc_func_list
        with open('/AlphaFold_refine/UltraDataset/Database1+/cc/cc_label.txt', 'r') as f:
            for line in f:
                line = line.strip('\n')
                cc_func.append(line)
        for i in cc_dict.keys():
            label = np.zeros(len(cc_func))
            for j in cc_dict[i]:
                if j in cc_func:
                    index = cc_func.index(j)
                    label[index] = 1
            cc_label_dict[i] = label
        return cc_label_dict, len(cc_func)

def Get_eachGOterms_List():
    mf_func_list = []
    bp_func_list = []
    cc_func_list = []
    with open('/AlphaFold_refine/UltraDataset/Database1+/mf/mf_label.txt', 'r') as f:
        for line in f:
            line = line.strip('\n')
            mf_func_list.append(line)

    with open('/AlphaFold_refine/UltraDataset/Database1+/bp/bp_label.txt', 'r') as f:
        for line in f:
            line = line.strip('\n')
            bp_func_list.append(line)

    with open('/AlphaFold_refine/UltraDataset/Database1+/cc/cc_label.txt', 'r') as f:
        for line in f:
            line = line.strip('\n')
            cc_func_list.append(line)
    return mf_func_list,bp_func_list,cc_func_list

def Generate_pdb_mapping_seq(IDlist,savepath):
    """
    用来生成pdb_mapping_seq 文件 ，该文件的每一行格式如下：
    >5UKL-G_P59768	SIAQARKLVEQLKMEANIDRIKVSKAAADLMAYCEAHAKEDPLLTPVPASENPFREKKFFCA
    - P59768: UPID
    - 5UKL: PDBID
    - -G: 指选中的蛋白质链为G （一个蛋白质可以有多条链）
    P59768 蛋白质在uniprot数据库中检索可得完整的氨基酸序列：
    MASNNTASIAQARKLVEQLKMEANIDRIKVSKAAADLMAYCEAHAKEDPLLTPVPASENPFREKKFFCAIL （由71个氨基酸残基组成）
    其在PDB数据库中对应的高精度空间结构：5UKL-G:
    - ID = 5UKL
    - method(使用的结构解析/测定技术) = X-ray ;
    - Resolution(分辨率)=2.15A ;
    - Chain(链) = G
    - Position = 8-69 (这说明这个蛋白质链是从完整序列的第8个氨基酸开始到第69个氨基酸)，
    即：MASNNTA-[SIAQARKLVEQLKMEANIDRIKVSKAAADLMAYCEAHAKEDPLLTPVPASENPFREKKFFCA]-IL
    括号里面的氨基酸残基序列才是5UKL-G对应的序列。
    :return:
    """
    # 第一步，提取蛋白质链对应position即[]中的序列信息
    startpoint = 0
    endpoint = 0
    pdb_mapping_seq_item_list = []
    # IDlist = ['2I0V$A=538-678, A=753-922$P07333','6KUW$A/B=29-241, A/B=372-458$P18825'] #用来debug
    for filteredID in IDlist: # 每一个filteredID的形式都为：1UMK$A=27-301$P00387，即以$分割
        # 使用 "$" 拆分字符串
        split_list = filteredID.split("$") # 此时split_list = ['1UMK', 'A=27-301', 'P00387']
        PDBID = split_list[0]
        UPID = split_list[2]
        # print(UPID)
        chaintype_and_range = split_list[1]  # 获取第二个元素，并进一步拆分
        if(len(chaintype_and_range.split('=')) > 2): #有的时候，链有不同的划分，如：2I0V$A=538-678, A=753-922$P07333 ，这个例子有两个”=“，得单独考虑
            chaintype_ = chaintype_and_range.split('=')[0][0]
            range_ = chaintype_and_range.split('=')[1].split(',')[0] #用'='划分后，继续用 ','划分第二个元素，划分后的第一个元素就是range(默认取第一个range)
            startpoint = int(range_.split('-')[0])  # 提取起始位置
            endpoint = int(range_.split('-')[1])  # 提取终止位置
        else:
            chaintype_ = chaintype_and_range.split('=')[0][0]  # 提取蛋白质链的类型,默认就选”A/B/C“中的第一个，即A
            range_ = chaintype_and_range.split('=')[1]  # 提取数值范围
            startpoint = int(range_.split('-')[0])  # 提取起始位置
            endpoint = int(range_.split('-')[1])  # 提取终止位置
        """根据UPID从ProGO中查询完整的氨基酸seq"""
        collection = MongoConnection()

        query_result = collection.find_one({"UPID": UPID})["Sequence"] # 从MonGoDB中通过UPID索引到其对应条目的 Sequence
        # print(UPID + "的IDENTIFIER_entryId 查询结果:", query_result, '\n返回类型：', type(query_result))
        seq = query_result
        chain_seq = seq[startpoint-1:endpoint] #截取从startpoint到endpoint的子串

        pdb_mapping_seq_item = ">" + str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID) + "\t" + str(chain_seq)
        pdb_mapping_seq_item_list.append(pdb_mapping_seq_item)
    List_to_TXT(testList=pdb_mapping_seq_item_list,savePath=savepath) #"G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/pdb_mapping_seq_new.txt"

def Generate_Go_label(IDlist,savepath):
    """
    用来生成 Go_label.txt,该文件的格式如下：
    以测试数据 即：UPID=P59768 为例：
    {'mf': ['GO:0070330', 'GO:0020037', 'GO:0005506', 'GO:0016712'], 'bp': ['GO:0019369', 'GO:0042178', 'GO:0006805'], 'cc': ['GO:0005737', 'GO:0043231', 'GO:0016020', 'GO:0005739']}

    于是它在Go_label.txt中：
    >5UKL-G_P59768（AF-A0A087X1C5-F1 这个ID在MongoDB中的Structure_AlphaFold条目下可以找到,可以通过这个ID爬取AlphaFold的PDB文件）
    mf:GO:0070330   GO:0020037  GO:0005506  GO:0016712
    bp:GO:0019369   GO:0042178  GO:0006805
    cc:GO:0005737   GO:0043231  GO:0016020  GO:0005739
    :param IDlist:
    :return:
    """
    with open(savepath, 'w') as f:
        for filteredID in IDlist:  # 每一个filteredID的形式都为：1UMK$A=27-301$P00387，即以$分割
            # 使用 "$" 拆分字符串
            split_list = filteredID.split("$")  # 此时split_list = ['1UMK', 'A=27-301', 'P00387']
            PDBID = split_list[0]
            UPID = split_list[2]
            chaintype_and_range = split_list[1]  # 获取第二个元素,即蛋白质链标识符以及position，并进一步拆分
            chaintype_ = chaintype_and_range.split('=')[0][0] #拿到蛋白质标识符
            intact_ID =  ">" + str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID)
            print("intact_ID:",intact_ID)
            """根据UPID从ProGO中查询完整的氨基酸 Gene ontology 标签"""
            # 先拿到原数据集中的mf,bp,cc标签，我们扩充的蛋白质数据的mf,bp,cc标签需要是和源数据集的交集，这样可以确保扩充数据集的同时避免标签的增加！
            mf_label_list_raw,bp_label_list_raw,cc_label_list_raw = Get_eachGOterms_List()
            # print(len(mf_label_list_raw),mf_label_list_raw,'\n',
            #       len(bp_label_list_raw),bp_label_list_raw,'\n',
            #       len(cc_label_list_raw),cc_label_list_raw)
            collection = MongoConnection(collection_name = 'ProGO_Target_Info')
            mf_list = collection.find_one({"UPID": UPID})["mf"]  # 从MonGoDB中通过UPID索引到其对应条目的 Golabel
            bp_list = collection.find_one({"UPID": UPID})["bp"]
            cc_list = collection.find_one({"UPID": UPID})["cc"]
            # 求与原标签列表的交集
            common_mf = set(mf_list) & set(mf_label_list_raw)
            common_bp = set(bp_list) & set(bp_label_list_raw)
            common_cc = set(cc_list) & set(cc_label_list_raw)
            common_mf = list(common_mf)
            common_bp = list(common_bp)
            common_cc = list(common_cc)
            print(len(common_mf), common_mf, type(common_mf),'\n',
                  len(common_bp), common_bp, type(common_bp),'\n',
                  len(common_cc), common_cc, type(common_cc))
            """开始写入 Go_label.txt文件
            每一轮for循环都在文件中写入一个intactID对应的GO标签，如：
            >1UMK-A_P00387
            mf: #列表1中的元素，以tab分割（列表为空就只保留一个mf:）
            bp: #列表2中的元素，以tab分割
            cc: #列表3中的元素，以tab分割
            """
            # with open(savepath, 'w') as f:
            # 先写入intactID
            f.write(f"{intact_ID}\n")
            # 写入mf交集
            mf_line = "mf:\t" + "\t".join(common_mf) if common_mf else "mf:"
            f.write(f"{mf_line}\t\n")

            # 写入bp交集
            bp_line = "bp:\t" + "\t".join(common_bp) if common_bp else "bp:"
            f.write(f"{bp_line}\t\n")

            # 写入cc交集
            cc_line = "cc:\t" + "\t".join(common_cc) if common_cc else "cc:"
            f.write(f"{cc_line}\t\n")

        # break #测试

def PDBFile_scrapy(IDlist,savepath,Structural_source='PDB'):
    """根据指定的PDBID,UPID,AlphaFold_indexID 爬取对应的PDB结构文件"""

    if not os.path.exists(savepath):
        os.makedirs(savepath)
    # test = 0

    count = 0
    for filteredID in IDlist:
        # 使用 "$" 拆分字符串
        split_list = filteredID.split("$")  # 此时split_list = ['1UMK', 'A=27-301', 'P00387']
        PDBID = split_list[0]
        if (PDBID in ["5VYC"]):
            print(f"该{PDBID}蛋白质结构属于大分子结构，只能下载为PDBx/mmCIF 文件")
            continue
        UPID = split_list[2]
        chaintype_and_range = split_list[1]  # 获取第二个元素，并进一步拆分
        chaintype_ = chaintype_and_range.split('=')[0][0]  # 提取蛋白质链的类型,默认就选”A/B/C“中的第一个，即A

        if (len(chaintype_and_range.split(
                '=')) > 2):  # 有的时候，链有不同的划分，如：2I0V$A=538-678, A=753-922$P07333 ，这个例子有两个”=“，得单独考虑
            chaintype_ = chaintype_and_range.split('=')[0][0]
            range_ = chaintype_and_range.split('=')[1].split(',')[
                0]  # 用'='划分后，继续用 ','划分第二个元素，划分后的第一个元素就是range(默认取第一个range)
            startpoint = int(range_.split('-')[0])  # 提取起始位置
            endpoint = int(range_.split('-')[1])  # 提取终止位置
        else:
            chaintype_ = chaintype_and_range.split('=')[0][0]  # 提取蛋白质链的类型,默认就选”A/B/C“中的第一个，即A
            range_ = chaintype_and_range.split('=')[1]  # 提取数值范围
            startpoint = int(range_.split('-')[0])  # 提取起始位置
            endpoint = int(range_.split('-')[1])  # 提取终止位置


        if (Structural_source == 'PDB'):
            pdb_scrapy_url =  "https://files.rcsb.org/download/"+ str(PDBID) + ".pdb" # 接口链接：https://files.rcsb.org/download/4WHJ.pdb  这里的4WHJ就是PDBID
            scrapy_url = pdb_scrapy_url
            File_name = str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID)
            # Filename这样设定，为了方便后续将这个.pdb文件转换为同名的.npy文件（接触图）即：5UKL-G_P59768_startpoint_endpoint.pdb -> 5UKL-G_P59768.npy
        else: #说明是要爬AlphaFold里面的预测蛋白质结构了,AlphaFold 数据查询接口：https://alphafold.ebi.ac.uk/files/"UPID"
            """AlphaFold的最新结构数据下载链接去ProGO中查询"""
            collection = MongoConnection()
            query_result = collection.find_one({"UPID": UPID})["Structure_AlphaFold"]["PDBfileDownloadUrl"]
            alphafold_scrapy_url = str(query_result)
            scrapy_url = alphafold_scrapy_url
            File_name = str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID) + "_" + str(
                startpoint) + "_" + str(endpoint) #此时的filename h爱需要再加上start point 和 end point

        try:
            # 发送请求并下载文件
            response = requests.get(scrapy_url)
            response.raise_for_status()  # 检查请求是否成功

            # 构建保存文件的路径，文件名以关键字命名
            file_path = os.path.join(savepath, f'{File_name}.pdb')  #这里的 .pdb 可以根据文件类型修改

            # 写入文件
            with open(file_path, 'wb') as file:
                file.write(response.content)
            count += 1
            print(f"文件{count}： {file_path} 下载成功")
            # if(test == 1):
            #     print("测试结束")
            #     break

        except requests.exceptions.RequestException as e:
            print(f"下载Structural_source={Structural_source}=={scrapy_url} 时出错: {e}")
            break

# 定义tensorflow解析函数
def parse_function(example_proto):
    # 定义特征描述，以便解析
    feature_description = {
        'feature1': tf.io.FixedLenFeature([], tf.float32),
        'label': tf.io.FixedLenFeature([], tf.int64),
    }
    # 解析样本
    return tf.io.parse_single_example(example_proto, feature_description)

def tensorflowDataReader(tfrecord_file):
    """定义解析函数，用于解析DeepFRI提供的TFRecord数据集文件中的每个样本"""
    # 使用TFRecordDataset读取文件
    raw_dataset = tf.data.TFRecordDataset(tfrecord_file)
    parsed_dataset = raw_dataset.map(parse_function) # 解析数据集
    # 查看解析后的内容
    for parsed_record in parsed_dataset:
        print("Feature1:", parsed_record['feature1'].numpy())
        print("Label:", parsed_record['label'].numpy())


def statistics_ProGO(ProGO_version = '3.0',testcount = 47789):
    """
    统计ProGO数据库
    :return:
    """
    collection = MongoConnection(collection_name='ProGO_Target_Info')
    statistics_collection = MongoConnection(collection_name='ProGO_statistics_Info')
    cursor = collection.find()  # 获取一个游标（cursor），然后逐个迭代游标以按顺序遍历MongoDB中的每个条目访问每个文档。
    record = {}
    PDB_statistics_list = []
    PDB_statistics_info = {}
    without_alpha_list = []
    MFlabel_statistics_list = []
    BPlabel_statistics_list = []
    CClabel_statistics_list = []
    MFlabel_statistics_info = {}
    BPlabel_statistics_info = {}
    CClabel_statistics_info = {}
    inactive_count = 0
    errorID = ["Q54436", "Q6FNA9", "I6Y8B5", "A0A1C3NEV1", "A5K421", "Q8GRC7",
               "Q7RSE5", "M7TVE7", "Q64823", "Q30JB5",
               "Q6KC22"]  # 这些ID是在爬虫阶段产生错误的ID无法解析他们的数据接口返回值，需要后续单独用beautifulsoup或其他方法直接解析他们的网页HTML拿数据
    protein_entries_count = 0

    for document in cursor:
        if (testcount == protein_entries_count):
            break
        protein_entries_count += 1
        UPID = document["UPID"]
        print(f"当前检索条目为:{protein_entries_count}","UPID=",UPID)
        if (document["UPID"] in errorID):
            print("错误ID", document["UPID"])
            continue
        if ('inactive' in document):
            # 数据库中还存在 未激活状态的数据比如：
            # {'_id': ObjectId('662deb4021712c0f197615b9'), 'UPID': 'Q08AF8', 'inactive': 'True', 'inactiveReasonType': ['DEMERGED'], 'bp': [], 'cc': [], 'mf': []}
            print("该数据在UniprotKB中未激活，inactiveReasonType = ", document['inactiveReasonType'])
            inactive_count += 1
            continue
        if (getAlphaEntryIDFromMongo(document["UPID"]) == 'null'):  # 在后续的测试中发现，有些蛋白质数据是没有AlphaFold预测结构的（比如UPID =A2VEC9、P22352..... ），这会导致数据缩水，所以要提前排除
            print("该UPID=", document["UPID"], "的蛋白质没有与之对应的AlphaFold结构")
            without_alpha_list.append(document["UPID"])
            continue
        if (len(document["Structure_PDB"]) == 0):  # 筛选出有高质量空间结构即PDBID的蛋白质条目
            continue

        for PDB_structure_json in document["Structure_PDB"]:  # 遍历所有PDB结构，选出我们要的
            PDBID = PDB_structure_json['id']
            if PDBID not in PDB_statistics_list:
                PDB_statistics_list.append(PDBID)

        # try:
        #     testlist = len(document["mf"])
        #     testlist = len(document["bp"])
        #     testlist = len(document["cc"])
        # except:
        #     print(f"该条目{document}没有添加单独的mf,bp,cc信息")
        #     continue

        if (len(document["mf"]) != 0): #这里的document["mf"] 可能为空列表
            for mfterm in document["mf"]:
                if mfterm not in MFlabel_statistics_list:
                    MFlabel_statistics_list.append(mfterm)
        if (len(document["bp"]) != 0): #这里的document["bp"] 可能为空列表
            for bpterm in document["bp"]:
                if bpterm not in BPlabel_statistics_list:
                    BPlabel_statistics_list.append(bpterm)
        if (len(document["cc"]) != 0): #这里的document["cc"] 可能为空列表
            for ccterm in document["cc"]:
                if ccterm not in CClabel_statistics_list:
                    CClabel_statistics_list.append(ccterm)
    MFlabel_statistics_info['count'] = len(MFlabel_statistics_list)
    MFlabel_statistics_info['info'] = MFlabel_statistics_list

    BPlabel_statistics_info['count'] = len(BPlabel_statistics_list)
    BPlabel_statistics_info['info'] = BPlabel_statistics_list

    CClabel_statistics_info['count'] = len(CClabel_statistics_list)
    CClabel_statistics_info['info'] = CClabel_statistics_list

    PDB_statistics_info['count'] = len(PDB_statistics_list)
    PDB_statistics_info['info'] = PDB_statistics_list



    record['ProGO_versionID'] = ProGO_version
    record['PDB_entries_statistics'] = PDB_statistics_info
    record['without_alpha_entries'] = without_alpha_list
    record['MFlabel_statistics'] = MFlabel_statistics_info
    record['BPlabel_statistics'] = BPlabel_statistics_info
    record['CClabel_statistics'] = CClabel_statistics_info
    record['errorJson_ID'] = errorID
    record['protein_entries_count'] = protein_entries_count
    record['inactive_entries_count'] = inactive_count
    record['statistics_num'] = testcount

    statistics_collection.update_one({'ProGO_versionID': record['ProGO_versionID'],'statistics_num': record['statistics_num']}, {'$set': record},
                                     True)  # 存储没有的ProGO版本对应条目，或更新现有的



if __name__ == '__main__':
    """preprocessing"""
    # Gen_GO_labelTXT() #给Mongo中的数据生成mf,bp,cc索引，起初只给Mongo数据库中前5000条数据做了这个更新
    # LabelhaveMf_andAll_InMongo() #返回筛选出来的ID列表,并将其存为同级目录下的的intact_label.txt
    # GenIntact_labelTXT() #生成最终筛选的 AlphaFoldentryID与UPID的组合的Intact_label.txt

    """使用LabelhaveMf_andAll_InProGO()从ProGO中筛选除符合条件的UPID"""
    # mf_label_list,bp_label_list,cc_label_list = Get_eachGOterms_List()
    # print(len(mf_label_list),mf_label_list,'\n',
    #       len(bp_label_list),bp_label_list,'\n',
    #       len(cc_label_list),cc_label_list)
    #
    # IntactLabelIDlist = LabelhaveMf_andAll_InProGO(filterNum=302,mf_func_list=mf_label_list,bp_func_list=bp_label_list,cc_func_list=cc_label_list)
    # print("筛选出符合条件的ID列表：",IntactLabelIDlist)
    # List_to_TXT(testList=IntactLabelIDlist,savePath="G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")

    """基于filteredIDlist生成扩展的pdb_mapping_seq_new.txt"""
    # filteredIDlist = TXT_to_List(filePath = "G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")
    # Generate_pdb_mapping_seq(IDlist = filteredIDlist,savepath= "G:/fkfk/AlphaFold_refine/UltraDataset/testData3_290entries/pdb_mapping_seq_new.txt")
    # #savepath = G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/pdb_mapping_seq_new.txt

    """基于filteredIDlist生成扩展的Go_label_new.txt"""
    # filteredIDlist = TXT_to_List(filePath="G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")
    # Generate_Go_label(IDlist=filteredIDlist,savepath = "G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/Go_label_new.txt")

    """基于filteredIDlist 中的PDBID爬取(下载对应的PDB文件)到指定路径"""
    # filteredIDlist = TXT_to_List(filePath="G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")
    # PDBFile_scrapy(IDlist=filteredIDlist,savepath = "G:/fkfk/AlphaFold_refine/biotoolbox/PDBfile_Dataset_Ultimate",
    #                Structural_source='PDB') #savepath = G:\fkfk\AlphaFold_refine\biotoolbox\PDBfile_Dataset_Ultimate

    """基于filteredIDlist 中的UPID去MongoDB中索引AlphaFold_index,之后再AlphaFold爬虫(下载对应的PDB文件)到指定路径"""
    # # 注意!若此处发生网络接口错误，记得开VPN
    # filteredIDlist = TXT_to_List(filePath="G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")
    # PDBFile_scrapy(IDlist=filteredIDlist,
    #                savepath="G:/fkfk/AlphaFold_refine/biotoolbox/AlphaFoldfile_Dataset_Ultimate",
    #                Structural_source='AlphaFold')  # PDBfile 的 savepath = G:\fkfk\AlphaFold_refine\biotoolbox\PDBfile_Dataset_Ultimate

    """生成290条 train_data_extend(即mf,bp,cc文件夹下的train_data 文件)"""
    # filteredIDlist = TXT_to_List(filePath="G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")
    # with open("G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/mf/train_data_new.txt", 'w') as file:
    #     for filteredID in filteredIDlist:  # 每一个filteredID的形式都为：1UMK$A=27-301$P00387，即以$分割
    #         # 使用 "$" 拆分字符串
    #         split_list = filteredID.split("$")  # 此时split_list = ['1UMK', 'A=27-301', 'P00387']
    #         PDBID = split_list[0]
    #         UPID = split_list[2]
    #         chaintype_and_range = split_list[1]  # 获取第二个元素,即蛋白质链标识符以及position，并进一步拆分
    #         chaintype_ = chaintype_and_range.split('=')[0][0]  # 拿到蛋白质标识符
    #         intact_ID = str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID)
    #         file.write(intact_ID + "\n")


    """接下来，最后一步，通过上面爬取的PDB文件，生成contact map"""
    # contact_mapping_PDB(PDBID='3UDW-C_P15151', threshold=8., GCN_transformer=1,savepath="G:/fkfk/AlphaFold_refine/biotoolbox/Contactmapfile_Dataset_Ultimate/")
    # filteredIDlist = TXT_to_List(filePath="G:/fkfk/AlphaFold_refine/UltraDataset/Dataset_Ultimate/filteredIDlist.txt")
    # count = 0
    # error_id_list = []  # 存放 下载来的PDB文件有格式错误的PDBID
    # for filteredID in filteredIDlist:  # 每一个filteredID的形式都为：1UMK$A=27-301$P00387，即以$分割
    #     # 使用 "$" 拆分字符串
    #     split_list = filteredID.split("$")  # 此时split_list = ['1UMK', 'A=27-301', 'P00387']
    #     PDBID = split_list[0]
    #     UPID = split_list[2]
    #     chaintype_and_range = split_list[1]  # 获取第二个元素,即蛋白质链标识符以及position，并进一步拆分
    #     chaintype_ = chaintype_and_range.split('=')[0][0]  # 拿到蛋白质标识符
    #     intact_ID = str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID)
    #
    #     try:
    #         contact_mapping_PDB(PDBID=str(intact_ID),threshold=8.,GCN_transformer=1,savepath = "G:/fkfk/AlphaFold_refine/biotoolbox/Contactmapfile_Dataset_Ultimate/")
    #         count += 1
    #         print(f"intact_ID={intact_ID}的蛋白质条目接触图已经生成，{count}")
    #     except:
    #         print(f"filteredID={PDBID},条目的PDB文件格式错误，下一个")
    #         error_id_list.append(PDBID)
    #         continue
    # print(error_id_list)
    # pass

    """分析PDB-GO（DeepFRI提供的tensorflow 数据集）"""
    # tensorflowDataReader('G:/PDB-GO/PDB_GO_train_01-of-30.tfrecords')

    """对完整的ProGO数据库进行统计分析"""
    ProGOVersion = input("输入需统计的ProGO版本（默认3.0）")
    Testnumber = int(input("输入您需要统计ProGO数据库的前n个条目（默认为统计完整ProGO数据库,i.e ProGO 3.0 entries=47789）"))
    statistics_ProGO(ProGO_version=ProGOVersion,testcount=Testnumber)







