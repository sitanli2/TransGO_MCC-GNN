import CollectFromDisk as CFD
from AlphaFold_refine.Protein_spider.UniProt_Scrapy.uniprot_refineFinal import mainscrapy

if __name__ == '__main__':
    # """得到要爬取的UPID列表 并村存入txt文件中去"""
    # FullUPIDlist,FullPIDlist = CFD.GetFullIDlist()
    # with open('IDlist.txt','w') as file:
    #     for item in FullUPIDlist:
    #         file.write(item + '\n')

    UPIDlist = []
    with open('IDlist.txt', 'r') as file:
        for line in file:
            UPIDlist.append(line.strip()) # 去掉每行末尾的换行符，并添加到新的列表中
    #print("UPID=",UPIDlist[2518])
    mainscrapy(UPIDlist)
