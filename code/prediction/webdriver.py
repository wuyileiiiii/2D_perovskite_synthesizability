#Based on the CID of molecules, this code can automatically download 3D structures of molecules from PubChem
#Require file "CID.txt" with the CID of molecules
#Require file "fake_useragent_0.1.11.json"
#The inconsistence of the version between Webdriver and IE might cause errors
#Solutions: Find the version of IE in folder "C:\Program Files (x86)\Microsoft\Edge\Application"/
#Then download the corresponding version of webdriver in "https://developer.microsoft.com/en-us/microsoft-edge/tools/webdriver/",and change file "msedgewebdriver.exe"
#author: ylwu
#Date: 2022/10/7



import re
import os
from fake_useragent import UserAgent
from selenium import webdriver
import time

browser = webdriver.Edge('C:\Program Files (x86)\Microsoft\Edge\Application\msedgedriver.exe')
url='https://pubchem.ncbi.nlm.nih.gov/compound/cid'
CID=[]
location = os.getcwd() + '/fake_useragent_0.1.11.json'
ua=UserAgent(path=location)
CID_file=open("CID.txt",'r',encoding='utf-8')
cid_line=CID_file.readlines()

def ttt(url):
    browser.get(url)

print(len(cid_line))
url1='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/cid/record/SDF/?record_type=3d&amp;response_type=save&amp;response_basename=Conformer3D_CID_cid'
for i in range(len(cid_line)):
    cid_url1=re.sub(r'cid',cid_line[i],url1)
   # Replace the CAS in the url
    headers = {"User-Agent": ua.random}
    print(cid_url1)
    ttt(cid_url1)
    text=browser.page_source

browser.close()