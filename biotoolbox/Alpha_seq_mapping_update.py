"""
last scrapy loop output:
intactID = 5ZJR-A_P09087, start point = 369,end point=452,location = A=369-452, AlphaFileName = 5ZJR-A_P09087_369_452
第7177个文件已成功保存到: G:/fkfk/AlphaFold_refine/biotoolbox/AlphaFoldfile_Dataset_Ultimate/5ZJR-A_P09087_369_452.pdb

"""
import os
from pymongo import MongoClient



def PDBFile_bug_detection():
    seq_error_list = ['1G91-A_P55773', '6C95-D_Q9NX55', '2GK1-I_P06865', '1LO1-A_O95718', '5M73-D_O76094', '5G04-X_Q9UJX3', '6FOH-A_Q6P587', '2V6C-A_P50580', '3ABR-B_P19636', '1E5W-A_P26038', '4UHT-A_P0AE88', '2OOY-G_Q10343', '4B93-A_P70280', '3EZZ-A_Q13115', '1QGP-A_P55265', '2LEF-A_P27782', '5HGZ-A_Q9H7X0', '2MF8-A_Q8CFC2', '1F5X-A_P27870', '2MGY-A_P50637', '2F5Z-K_O00330', '2L26-A_P9WIU5', '2YF0-A_Q9Y217', '2WYA-A_P54868', '4RLT-A_P9WFK1', '1G5V-A_Q16637', '1QU9-A_P0AF93', '2GK1-N_P07686', '2PFU-A_P0ABV2', '4G8X-A_P0A0J0', '3KYD-D_P63165', '3GBJ-A_Q9NQT8', '1TWI-A_Q58497', '1DI7-A_P0AF03', '5DAR-B_P54049', '5V4P-A_P41816', '3AUW-B_P48542', '2VM5-A_Q13075', '2FSE-A_P01903', '1EQT-A_P13501', '3ZZW-A_Q01974', '2EAX-A_Q96LB8', '1BFV-L_P01631', '2WO3-B_O43921', '5HYD-A_Q8WXG8', '1JPE-A_P36655', '2PG7-A_Q16696', '1W2F-A_P23677', '5DJO-A_Q9EQW7', '2M3D-A_Q9NR30', '2MZW-A_Q2G0N1', '6C9M-B_P41227', '1E41-A_Q13158', '1A0A-A_P07270', '2KD3-A_Q99P68', '1LV7-A_P0AAI3', '1TUZ-A_P23743', '1DCF-A_P49333', '5T7H-A_P00044', '1F9A-A_Q57961', '4XT3-B_P78423', '1SXL-A_P19339', '1JBI-A_O43405', '1W3B-A_O15294', '4JPZ-A_Q92913', '2D68-A_O95684', '2V0V-A_Q14995', '2H6U-A_Q06S87', '2L3L-A_Q15814', '5KZQ-A_Q14416', '1O6X-A_P48052', '6H0G-B_Q96SW2', '2ZJ1-A_P9WGV3', '1IUY-A_Q9JLV5', '5U1S-B_P40316', '1DQB-A_P07204', '3DXS-X_Q9S7J8', '4BHB-A_P9WJW5', '4PD4-K_P01647', '1XWB-A_Q9V429', '1PQV-D_P20433', '1XJH-A_P0A6Y5', '6CQ6-A_P97438', '2NLL-B_P10828', '1JLI-A_P08700', '3PYF-A_P9WKG7', '1G8Q-A_P60033', '4RH8-A_P0ADC1', '2L3X-A_Q6H8Q1', '2EW9-A_P35670', '1VND-A_P22808', '1BR0-A_P32851', '1GPQ-A_P0AD59', '2WFD-A_Q9P2J5']#对应的Alpha_seq_mapping 序列做修正
    error_dict = {'1G91-A_P55773': 'Seq_mapping_error', '6C95-D_Q9NX55': 'Seq_mapping_error', '2GK1-I_P06865': 'Seq_mapping_error', '4ML8-A_Q709Q5': ['TypeError', TypeError('string indices must be integers')], '1LO1-A_O95718': 'Seq_mapping_error', '4N6B-A_I1KHY6': ['TypeError', TypeError('string indices must be integers')], '5M73-D_O76094': 'Seq_mapping_error', '5G04-X_Q9UJX3': 'Seq_mapping_error', '6FOH-A_Q6P587': 'Seq_mapping_error', '2V6C-A_P50580': 'Seq_mapping_error', '3ABR-B_P19636': 'Seq_mapping_error', '1E5W-A_P26038': 'Seq_mapping_error', '4UHT-A_P0AE88': 'Seq_mapping_error', '4BUM-X_Q8AWD0': ['TypeError', TypeError('string indices must be integers')], '2OOY-G_Q10343': 'Seq_mapping_error', '4B93-A_P70280': 'Seq_mapping_error', '3EZZ-A_Q13115': 'Seq_mapping_error', '1QGP-A_P55265': 'Seq_mapping_error', '2LEF-A_P27782': 'Seq_mapping_error', '5HGZ-A_Q9H7X0': 'Seq_mapping_error', '2MF8-A_Q8CFC2': 'Seq_mapping_error', '1F5X-A_P27870': 'Seq_mapping_error', '2MGY-A_P50637': 'Seq_mapping_error', '4AGS-A_A4I8P2': ['TypeError', TypeError('string indices must be integers')], '2F5Z-K_O00330': 'Seq_mapping_error', '2L26-A_P9WIU5': 'Seq_mapping_error', '2YF0-A_Q9Y217': 'Seq_mapping_error', '2WYA-A_P54868': 'Seq_mapping_error', '4RLT-A_P9WFK1': 'Seq_mapping_error', '4L8Y-C_Q6NY42': ['TypeError', TypeError('string indices must be integers')], '1G5V-A_Q16637': 'Seq_mapping_error', '6D21-A_E7FE33': ['TypeError', TypeError('string indices must be integers')], '1QU9-A_P0AF93': 'Seq_mapping_error', '2GK1-N_P07686': 'Seq_mapping_error', '2PFU-A_P0ABV2': 'Seq_mapping_error', '4G8X-A_P0A0J0': 'Seq_mapping_error', '3KYD-D_P63165': 'Seq_mapping_error', '1WN0-A_Q9SLX1': ['TypeError', TypeError('string indices must be integers')], '5JH1-A_A2T1W7': ['TypeError', TypeError('string indices must be integers')], '3GBJ-A_Q9NQT8': 'Seq_mapping_error', '1TWI-A_Q58497': 'Seq_mapping_error', '1DI7-A_P0AF03': 'Seq_mapping_error', '5DAR-B_P54049': 'Seq_mapping_error', '5V4P-A_P41816': 'Seq_mapping_error', '3AUW-B_P48542': 'Seq_mapping_error', '4LUR-A_F1Q9N9': ['TypeError', TypeError('string indices must be integers')], '2VM5-A_Q13075': 'Seq_mapping_error', '2M3C-A_Q5XTN3': ['TypeError', TypeError('string indices must be integers')], '2FSE-A_P01903': 'Seq_mapping_error', '1EQT-A_P13501': 'Seq_mapping_error', '3ZZW-A_Q01974': 'Seq_mapping_error', '5ZJI-G_B6U534': ['TypeError', TypeError('string indices must be integers')], '2EAX-A_Q96LB8': 'Seq_mapping_error', '1BFV-L_P01631': 'Seq_mapping_error', '2WO3-B_O43921': 'Seq_mapping_error', '5ZJI-K_B6TR16': ['TypeError', TypeError('string indices must be integers')], '5HYD-A_Q8WXG8': 'Seq_mapping_error', '1JPE-A_P36655': 'Seq_mapping_error', '2PG7-A_Q16696': 'Seq_mapping_error', '1W2F-A_P23677': 'Seq_mapping_error', '5DJO-A_Q9EQW7': 'Seq_mapping_error', '2M3D-A_Q9NR30': 'Seq_mapping_error', '2MZW-A_Q2G0N1': 'Seq_mapping_error', '6C9M-B_P41227': 'Seq_mapping_error', '1E41-A_Q13158': 'Seq_mapping_error', '1A0A-A_P07270': 'Seq_mapping_error', '2KD3-A_Q99P68': 'Seq_mapping_error', '1LV7-A_P0AAI3': 'Seq_mapping_error', '2KSZ-A_Q39890': ['TypeError', TypeError('string indices must be integers')], '1TUZ-A_P23743': 'Seq_mapping_error', '1DCF-A_P49333': 'Seq_mapping_error', '5T7H-A_P00044': 'Seq_mapping_error', '1F9A-A_Q57961': 'Seq_mapping_error', '4XT3-B_P78423': 'Seq_mapping_error', '1SXL-A_P19339': 'Seq_mapping_error', '1JBI-A_O43405': 'Seq_mapping_error', '6FKK-A_Q7KK54': ['TypeError', TypeError('string indices must be integers')], '1K8I-B_Q31099': ['TypeError', TypeError('string indices must be integers')], '1W3B-A_O15294': 'Seq_mapping_error', '4JPZ-A_Q92913': 'Seq_mapping_error', '2D68-A_O95684': 'Seq_mapping_error', '2V0V-A_Q14995': 'Seq_mapping_error', '2H6U-A_Q06S87': 'Seq_mapping_error', '2L3L-A_Q15814': 'Seq_mapping_error', '6HRV-A_Q6P4V5': ['TypeError', TypeError('string indices must be integers')], '5KZQ-A_Q14416': 'Seq_mapping_error', '1O6X-A_P48052': 'Seq_mapping_error', '6H0G-B_Q96SW2': 'Seq_mapping_error', '2ZJ1-A_P9WGV3': 'Seq_mapping_error', '1IUY-A_Q9JLV5': 'Seq_mapping_error', '5U1S-B_P40316': 'Seq_mapping_error', '1DQB-A_P07204': 'Seq_mapping_error', '2J3R-B_Q6DGL5': ['TypeError', TypeError('string indices must be integers')], '3DXS-X_Q9S7J8': 'Seq_mapping_error', '6B82-A_A2ATX9': ['TypeError', TypeError('string indices must be integers')], '4BHB-A_P9WJW5': 'Seq_mapping_error', '4PD4-K_P01647': 'Seq_mapping_error', '1XWB-A_Q9V429': 'Seq_mapping_error', '1PQV-D_P20433': 'Seq_mapping_error', '1XJH-A_P0A6Y5': 'Seq_mapping_error', '6CQ6-A_P97438': 'Seq_mapping_error', '6E0F-J_Q95U89': ['TypeError', TypeError('string indices must be integers')], '2NLL-B_P10828': 'Seq_mapping_error', '1JLI-A_P08700': 'Seq_mapping_error', '3PYF-A_P9WKG7': 'Seq_mapping_error', '1G8Q-A_P60033': 'Seq_mapping_error', '4RH8-A_P0ADC1': 'Seq_mapping_error', '2L3X-A_Q6H8Q1': 'Seq_mapping_error', '2EW9-A_P35670': 'Seq_mapping_error', '1VND-A_P22808': 'Seq_mapping_error', '1BR0-A_P32851': 'Seq_mapping_error', '2AQX-A_B2RXC2': ['TypeError', TypeError('string indices must be integers')], '6EPE-S_Q5U2S7': ['TypeError', TypeError('string indices must be integers')], '1GPQ-A_P0AD59': 'Seq_mapping_error', '2WFD-A_Q9P2J5': 'Seq_mapping_error'}
    print(len(seq_error_list))

    keys_with_seq_error = [key for key, value in error_dict.items() if value == 'Seq_mapping_error']
    print(len(keys_with_seq_error))

    keys_with_other_error = [key for key, value in error_dict.items() if (value != 'Seq_mapping_error')]

    print(len(keys_with_other_error),keys_with_other_error)

    for error_id in keys_with_other_error:
        print(error_id,':',error_dict[error_id])


    # 指定文件夹路径
    folder_path = 'G:/fkfk/AlphaFold_refine/biotoolbox/AlphaFoldfile_Dataset_Ultimate'  # 替换成你要查询的文件夹路径

    # 获取文件夹中所有文件的数量
    file_count = sum(1 for item in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, item)))
    file_name_list = [item for item in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, item))]


    # 打印文件个数
    # print(f"文件夹 '{folder_path}' 中的文件数量是: {file_count}")
    # print(f"文件夹 '{folder_path}' 中的文件列表: {file_name_list} 表长：{len(file_name_list)}")

    new_file_name_list = []

    for intactid in file_name_list:
        ID = "_".join(intactid.split('_')[:2])
        new_file_name_list.append(ID)

    # print(f"文件夹 '{folder_path}' 中的文件列表: {new_file_name_list} 表长：{len(new_file_name_list)}")


    Intactid_list_PDB_B = []
    dict_seq_len = {}
    with open('G:/fkfk/AlphaFold_refine/biotoolbox/' + 'alpha_mapping_seq.txt', 'r') as f:
        for line in f:
            if '>' in line:
                pid_chain_uid = line[1:].strip('\n').split('\t')[0]
                seq = line[1:].strip('\n').split('\t')[1]
                seq_len = int(len(seq))
                Intactid_list_PDB_B.append(pid_chain_uid)
                dict_seq_len[pid_chain_uid] = seq_len
    # print('\n',Intactid_list_PDB_B,len(Intactid_list_PDB_B))

    common_elements = list(set(new_file_name_list) & set(Intactid_list_PDB_B))

    different_elements = list(set(new_file_name_list) ^ set(Intactid_list_PDB_B))

    print(common_elements,len(new_file_name_list),len(Intactid_list_PDB_B),len(common_elements),len(different_elements))
    print(different_elements)
    return file_name_list

def List_to_TXT(testList,savePath): #换行写入
    with open(savePath,'w') as file:
        for item in testList:
            file.write(item + "\n")

def TXT_to_List(filePath): #换行读取成列表
    with open(filePath,'r') as file:
        IDlist = [line.strip() for line in file]
    # print(len(IDlist),IDlist)
    return IDlist

def MongoConnection(collection_name = 'ProGO_Target_Info'):
    mongoClient = MongoClient('mongodb://localhost:27017/')
    database = mongoClient.get_database('dtmr_dev')

    # collection = dtmr_dev['uniprot_Target_Info_test']

    collection = database.get_collection(collection_name)
    return collection

def Generate_alpha_mapping_seq(IDlist,savepath):
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
    # IDlist = ['6QEL-A_P0ACB0_1_471.pdb', '6QEL-G_P0AEF0_1_245.pdb', '6QI8-A_Q9Y265_1_456.pdb'] #用来debug
    for filteredID in IDlist: # 每一个filteredID的形式都为：1UMK$A=27-301$P00387，即以$分割
        # 拆分字符串 6QEL-A_P0ACB0_1_471.pdb
        split_list = filteredID.split("_") # 此时split_list = ['6QEL-A', 'P0ACB0', '1' , '471.pdb']

        PDBID = split_list[0]
        UPID = split_list[1]

        startpoint = int(split_list[2])
        endpoint = int(split_list[3].split('.')[0]) # '471.pdb' 拿到471


        """根据UPID从ProGO中查询完整的氨基酸seq"""
        collection = MongoConnection()

        query_result = collection.find_one({"UPID": UPID})["Sequence"] # 从MonGoDB中通过UPID索引到其对应条目的 Sequence
        # print(UPID + "的IDENTIFIER_entryId 查询结果:", query_result, '\n返回类型：', type(query_result))
        seq = query_result
        chain_seq = seq[startpoint-1:endpoint] #截取从startpoint到endpoint的子串

        pdb_mapping_seq_item = ">" + str(PDBID) +  "_" + str(UPID) + "\t" + str(chain_seq)
        pdb_mapping_seq_item_list.append(pdb_mapping_seq_item)
    List_to_TXT(testList=pdb_mapping_seq_item_list,savePath=savepath)

if __name__ == '__main__':
    """生成更新修订后的 alpha_mapping_seq_new.txt"""
    # FIle_nameList=PDBFile_bug_detection()
    # print(FIle_nameList,len(FIle_nameList))
    # Generate_alpha_mapping_seq(IDlist=FIle_nameList,savepath="G:/fkfk/AlphaFold_refine/biotoolbox/alpha_mapping_seq_new.txt")

    # 指定文件夹路径
    folder_path2 = 'G:/fkfk/AlphaFold_refine/biotoolbox/Contactmapfile_Dataset_Ultimate/contact_map_dense_alpha_8.0'
    # 获取文件夹中所有文件的数量
    file_count = sum(1 for item in os.listdir(folder_path2) if os.path.isfile(os.path.join(folder_path2, item)))
    # 打印文件个数
    print(f"文件夹 '{folder_path2}' 中的文件数量是: {file_count}")







