import requests
import pandas as pd
import re
import tqdm
#import requests
import urllib.request as request
from bs4 import BeautifulSoup
from jsonpath import jsonpath
from tkinter import _flatten #将二维列表转成一维列表
import time
from AlphaFold_refine.Protein_spider.UniProt_Scrapy.uniprot_refineFinal import mainscrapy

def FromPDBID_to_UPID(PDBID='4GRI',API = 0):
    """
    从一个指定的 PDB ID（蛋白质数据银行的条目ID）获取其对应的 UniProt ID (UPID)
    使用 RCSB PDB API(0)和 UniProt 的 RESTful API(1) 查询 PDB ID 和 UniProt ID 之间的映射关系
    :return:
    """
    if (API == 0): # RCSB PDB API
        url = f"https://data.rcsb.org/rest/v1/core/entry/{PDBID}"
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()
            uniprot_ids = []

            # 在结构实体中查找UniProt ID
            for entity in data.get('rcsb_entry_container_identifiers', {}).get('entity_polymer_entity', []):
                for ref in entity.get('rcsb_polymer_entity_container_identifiers', {}).get(
                        'reference_sequence_identifiers', []):
                    if ref['database_name'] == 'UniProt':
                        uniprot_ids.append(ref['database_accession'])

            return uniprot_ids if uniprot_ids else "No UniProt ID found."
        else:
            return f"Error: Could not retrieve data for PDB ID {PDBID}"

    elif (API == 1): # RESTful API
        url = "https://rest.uniprot.org/idmapping/run"
        headers = {
            "Content-Type": "application/x-www-form-urlencoded"
        }
        data = {
            "from": "PDB_ID",
            "to": "UniProtKB",
            "ids": PDBID
        }

        # 提交映射请求
        response = requests.post(url, headers=headers, data=data)

        if response.status_code == 200:
            job_id = response.json()["jobId"]

            # 检查作业状态，获取结果
            result_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
            result_response = requests.get(result_url)

            if result_response.status_code == 200 and result_response.json()["results"]:
                # 提取UniProt ID
                return result_response.json()["results"][0]["to"]
            else:
                return "No UniProt ID found or mapping failed."
        else:
            return "Error: Could not map PDB ID to UniProt ID."

    elif(API == 2): #SIFTS Mapping 接口
        success = 1
        url = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/" + str(PDBID)
        while True:
            try:
                response = requests.get(url)
                break
            except Exception as e:
                # 如果发生异常，则打印异常信息，并等待一段时间后重新尝试请求
                print("网页请求失败:", e)
                print("等待5秒后将重新尝试请求...")
                time.sleep(5)  # 等待5秒后重新尝试请求

        UPID_list = []

        if response.status_code == 200:
            data = response.json()
            print(f"query result of {PDBID}'s UPID(json)=",data)
            list1 = jsonpath(data,'$..UniProt')
            list2 = (list(_flatten(list1)))
            print(f"query list2={list2}")
            for jsondict in list2: #可能一个PDBID(比如说'5MRE') 映射多个UPID
                UPID_list = list(jsondict.keys())
                print(f"{PDBID}'s UPID mapping results={UPID_list}")
                break #只需要循环一次，避免后面有第二个json字典
            return UPID_list,success

        else:
            success = 0
            print(f"Error: Could not retrieve data for PDB ID {PDBID}")
            return UPID_list,success




def PDBinfo_scrapy():
    length = len(pd.read_csv("PDB_IDS_FOR_FILES.csv").keyword.values)

    pdb_ids = pd.read_csv("PDB_IDS_FOR_FILES.csv").keyword.values[1:]
    count = 0
    for pdb_id in tqdm.tqdm(pdb_ids):
        count += 1
        record = {}  # 定义一个集合
        html_url = "https://www.rcsb.org/structure/" + pdb_id
        response = request.urlopen(html_url)
        html = response.read()
        print('html:', html)
        soup = BeautifulSoup(html, 'html.parser')  # beautifulsoup(解析内容，解析器)
        # beautifulsoup是一个解析器，可以特定的解析出内容，省去了我们编写正则表达式进行正则查询的麻烦。这里指定Beautiful的解析器为“html.parser”
        # 编号，例如PROT000001
        code = code + 1
        record['code'] = 'PROT' + '0' * (6 - len(str(code))) + str(code)

        record['PDBID'] = pdb_id

        LigandsMainTable = soup.find("table", id="LigandsMainTable")


def GetUPID_Swiss_Prot():
    # 设置API的URL，过滤条件使用 reviewed:true 来只获取 Swiss-Prot 中的所有UPID条目
    url = "https://rest.uniprot.org/uniprotkb/stream?format=list&query=reviewed:true"

    # 发送请求并获取所有 Swiss-Prot 的 UniProt IDs
    response = requests.get(url)

    # 检查请求状态
    if response.status_code == 200:
        # 将所有UPID写入文件保存
        with open("G:/fkfk/AlphaFold_refine/Protein_spider/PDB_Scrapy/swissprot_uniprot_ids.txt", "w") as f:
            f.write(response.text)
        print("所有 Swiss-Prot UPIDs 已下载到 swissprot_uniprot_ids.txt 文件中。")
    else:
        print(f"请求失败，状态码: {response.status_code}")

def Generate_DeepFRI_mapping(endpoint):
    # 读取输入文件，并写入新的格式到输出文件
    # G:\fkfk\AlphaFold_refine\Protein_spider\PDB_Scrapy\PDB_id_DeepFRI.txt
    with open("G:/fkfk/AlphaFold_refine/Protein_spider/PDB_Scrapy/PDB_id_DeepFRI.txt", "r") as infile, \
            open("G:/fkfk/AlphaFold_refine/Protein_spider/PDB_Scrapy/UNIPROT_id_DeepFRI.txt", "a") as outfile: #这里要用a,不能用w，否则会覆盖原有的文件
        count = 0
        for line in infile:
            if (count <= endpoint): #跳过断点前的所有元素
                count += 1
                continue
            chainID = line.strip()  # 去除行末尾的换行符得到：2ZPA-A
            PDBID = chainID.split('-')[0] #去除链号，得到：2ZPA
            UPID_list,success_code = FromPDBID_to_UPID(PDBID=str(PDBID), API=2)  # "5MRE"(一个PDBID 映射多个UPID)
            print(f"result={UPID_list}")
            elements = ", ".join(UPID_list)  # 将列表中的元素用逗号分隔
            outfile.write(f"{chainID}:   {elements}\n")
            count += 1
            print(f"completed number={count}")
            # if (count==10):
            #     break

def Genrate_DeepFRI_upidlist():
    """
    根据上面生成的映射文件生成upid列表
    :return:
    """
    with open("G:/fkfk/AlphaFold_refine/Protein_spider/PDB_Scrapy/UNIPROT_id_DeepFRI.txt", "r") as file:
        Full_UPID_list = []
        for line in file: #4ZV0-B Q9I739, Q9I740
            # 去除两端空白并按冒号分割
            parts = line.strip().split(":")
            idlist = []
            if len(parts) > 1:
                # 获取冒号后面的部分，并按逗号分割成ID列表
                idlist = [id.strip() for id in parts[1].split(",")]
                # print(idlist)  # 输出当前行的ID列表
            else:
                print(f"{line}对应id为空列表")
            for id in idlist:
                if (id not in Full_UPID_list):
                    Full_UPID_list.append(id)
        return Full_UPID_list

def List_to_TXT(testList,savePath): #换行写入
    with open(savePath,'w') as file:
        for item in testList:
            file.write(item + "\n")
def TXT_to_List(filePath): #换行读取成列表
    with open(filePath,'r') as file:
        IDlist = [line.strip() for line in file]
    # print(len(IDlist),IDlist)
    return IDlist

if __name__ == '__main__':
    """获取Swiss_Prot下的所有UPID"""
    # GetUPID_Swiss_Prot()

    """基于DeepFRI提供的PDBID,查询其对应的UPID
    API = 2 使用SIFTS Mapping 接口，根据提供的PDBID,映射其对应的UPID"""
    # FromPDBID_to_UPID(PDBID='5MRE',API=2) # "5MRE"(一个PDBID 映射多个UPID)

    """根据DeepFRI提供的36641（29902 train , 3416 test , 3323 valid）条PDBID,映射所有的UPID"""
    # Generate_DeepFRI_mapping(endpoint=36641) #生成映射文件UNIPROT_id_DeepFRI.txt,输入断点 如4S3O-B=3824
    # Full_UPID_list = Genrate_DeepFRI_upidlist() #生成UPID list
    # print(Full_UPID_list)
    # print(f"len of the UPID_list no repeat = {len(Full_UPID_list)}")
    # List_to_TXT(Full_UPID_list,'G:/fkfk/AlphaFold_refine/Protein_spider/UniProt_Scrapy/DeepFRI_mapping_UPID')

    """开始ProGO爬虫"""
    uniprotId_list = TXT_to_List('G:/fkfk/AlphaFold_refine/Protein_spider/UniProt_Scrapy/DeepFRI_mapping_UPID')
    mainscrapy(UPID_list=uniprotId_list)



