# script to submit R generated fasta files to the mutpred2 form automatically

from selenium import webdriver
import chromedriver_binary
import time

email = "example@email.com
text = ">HSV1_UL23 L50M R51W Y53C Y53D D55N G56S P57H H58R H58N G61A G61E K62N T63A T63I T65N Q67R L68R L72V S74L R75C D76N D77N I78F E83K E83N P84S P84L M85I R89Q R89H A93V A98T I100T T103P H105P R106H Q109K E111K A118V M121R M121K Q125H G129D M130I P131S\nMASYPCHQHASAFDQAARSRGHNNRRTALRPRRQQEATEVRPEQKMPTLLRVYIDGPHGMGKTTTTQLLVALGSRDDIVYVPEPMTYWRVLGASETIANIYTTQHRLDQGEISAGDAAVVMTSAQITMGMPYAVTDAVLAPHIGGEAGSSHAPPPALTLIFDRHPIAALLCYPAARYLMGSMTPQAVLAFVALIPPTLPGTNIVLGALPEDRHIDRLAKRQRPGERLDLAMLAAIRRVYGLLANTVRYLQCGGSWREDWGQLSGTAVPPQGAEPQSNAGPRPHIGDTLFTLFRAPELLAPNGDLYNVFAWALDVLAKRLRSMHVFILDYDQSPAGCRDALLQLTSGMVQTHVTTPGSIPTICDLARTFAREMGEAN"

web = webdriver.Chrome()  # Optional argument, if not specified will search path
web.get('http://mutpred.mutdb.org/index.html')

time.sleep(5)

email_input = web.find_element_by_xpath('/html/body/div[1]/div/div/div/div/div[2]/form/div[1]/input')
email_input.send_keys(email)

fasta = text
fasta_input = web.find_element_by_xpath('/html/body/div[1]/div/div/div/div/div[2]/form/div[2]/textarea')
fasta_input.send_keys(fasta)


submit = web.find_element_by_xpath('/html/body/div[1]/div/div/div/div/div[2]/form/button')
submit.click()
