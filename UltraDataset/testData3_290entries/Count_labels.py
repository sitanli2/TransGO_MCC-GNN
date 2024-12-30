


def LabelhaveMf_andAll(getnum):
    """
    此函数用来从data_collect 那个Go_label原件中筛选出所有有MF标签的IDlist 即 mfIDlist
    随后再筛选出 BP,MF，CC 三个标签都有的IDlist 即 IntactLabelIDlist
    :parameter:getnum IntactLabelIDlist需要几条 ，即拥有完整bp,cc,mf 的ID
    :return: mfIDlist IntactLabelIDlist
    """
    IntactLabelIDlist = [] #用来存放mf,bp,cc三个标签都有的ID
    mfIDlist = [] #用来存放有mf标签的这些ID
    with open('../biLSTM+GCN/data_collect/Go_label.txt', 'r') as f: #data_collect 下的Go_label 原件
        for line in f:
            if '>' in line: #说明这一行是ID 如 >4WKG-A_P77398，也作为新的一个开始需要把所有Flag置为0
                mfflag = 0
                bpflag = 0
                ccflag = 0
                pdb_chain_uid = line[1:].strip('\n')
            elif 'mf' in line: #当前行是mf标签行
                if 'mf:' == line.strip('\n'): #说明这一整行 只有一个 ‘mf:’ ,即mf标签为空
                    # NoMfIDlist.append(pdb_chain_uid)
                    pass #当然也可以向上面一样存一下这个ID
                else:
                    mfflag= 1 #将flag置为1
                    mfIDlist.append(pdb_chain_uid)
            elif 'bp' in line: #当前行是bp标签
                if 'bp:' != line.strip('\n'): #看这一行bp标签是否为空，不为空将对应的flag置为一
                    bpflag = 1
            elif 'cc' in line: #当前行是cc标签,经过几次循环后下一行有 '>' 这个时候就要判断一下几个flag了
                if 'cc:' != line.strip('\n'): #看这一行cc标签是否为空，不为空将对应的flag置为一
                    ccflag = 1
                if (mfflag == 1 and bpflag == 1 and ccflag == 1): # 三个标签都有
                    IntactLabelIDlist.append(pdb_chain_uid)
            if(len(IntactLabelIDlist) >= getnum):
                break
        return mfIDlist,IntactLabelIDlist



def getnorepeat_mffuclist():
    """
    此函数用来返回testData下面的 Golabel文件中所有的ID中所含的所有不重复的mf标签，
    然后再将其通过下面的函数 generate_mf_labelTXT 写入testData下面的/mf/mf_label.txt
    :return:
    """
    mf_func = []
    mf_func_norepeat = []
    with open('./Go_label.txt', 'r') as f: #testdata 下面的Golabel

        for line in f:
            if '>' in line:
                pdb_chain_uid = line[1:].strip('\n')
            elif 'mf' in line:
                if 'mf:' == line.strip('\n'): #mf 标签为空
                    continue
                else:
                    mf_func_list = line[3:].strip().split('\t')
                    #print(mf_func_list,"\n")
                    for i in mf_func_list:
                        if i not in mf_func_norepeat:
                            mf_func_norepeat.append(i)
    return mf_func_norepeat

def generate_mf_labelTXT(mf_func_norepeat_file):
    """
    写入 未重复的 mf_label
    :param mf_func_norepeat_file:
    :return:
    """
    with open('./mf/mf_label.txt', 'w') as file:
        for item in mf_func_norepeat_file:
            file.write("%s\n" % item)


def CopyContactMapFromFilePath(path,filename_list,destination_folder):
    """

    :param path: 就是放所有contact map 的绝对路径
    :param filenames: 需要复制到新文件夹（testData）下的contact map的索引（文件名）列表,比如说filenames = [1A0A-A_P07270.npy,.....]
    :parameter destination_folder: testData下面放contact map 的文件夹
    :return:
    """
    import os
    import shutil
    # 确保目标文件夹存在，如果不存在则创建
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    for filename in filename_list:
        source_path = os.path.join(path, filename)
        destination_path = os.path.join(destination_folder, filename)
        shutil.copy(source_path, destination_path)
        print(f"文件 {filename} 已复制到 {destination_folder}")

def TXT_to_List(filePath): #换行读取成列表
    with open(filePath,'r') as file:
        IDlist = [line.strip() for line in file]
    # print(len(IDlist),IDlist)
    return IDlist

def GetIntactID_file():
    filteredIDlist = TXT_to_List(filePath="G:/fkfk/AlphaFold_refine/UltraDataset/testData3_290entries/filteredIDlist.txt")
    with open("G:/fkfk/AlphaFold_refine/UltraDataset/testData3_290entries/mf/train_data.txt", 'w') as file:
        for filteredID in filteredIDlist:  # 每一个filteredID的形式都为：1UMK$A=27-301$P00387，即以$分割
            # 使用 "$" 拆分字符串
            split_list = filteredID.split("$")  # 此时split_list = ['1UMK', 'A=27-301', 'P00387']
            PDBID = split_list[0]
            UPID = split_list[2]
            chaintype_and_range = split_list[1]  # 获取第二个元素,即蛋白质链标识符以及position，并进一步拆分
            chaintype_ = chaintype_and_range.split('=')[0][0]  # 拿到蛋白质标识符
            intact_ID = str(PDBID) + "-" + str(chaintype_) + "_" + str(UPID)
            file.write(intact_ID + "\n")



if __name__ == '__main__':
    """生成不重复的mf_label文件"""
    mf_func_norepeat_file = getnorepeat_mffuclist()
    print("no repeat func list:", mf_func_norepeat_file)
    generate_mf_labelTXT(mf_func_norepeat_file)
    """生成train_data.txt"""
    # GetIntactID_file()

    # i = input("输入你要的拥有完整label 的ID条数=")
    # mfIDlist,IntactIdlist = LabelhaveMf_andAll(int(i))
    # print(mfIDlist,'\n',IntactIdlist)
    # print(len(mfIDlist),'\n',len(IntactIdlist))

