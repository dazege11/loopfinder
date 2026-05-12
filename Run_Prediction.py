import argparse
import pandas as pd
from Bio import pairwise2
from Bio.Align import substitution_matrices
import esm
import torch
import torch.nn as nn
import numpy as np
import pickle
import importlib
import __main__
from tqdm import tqdm
import time


def ProteinEmbedding(pbar,data):
    model, alphabet = esm.pretrained.load_model_and_alphabet("esm2_t33_650M_UR50D.pt")
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]
    sequence_representations = []
    for i, (_, seq) in enumerate(data):
        sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0))
    pbar.update(3)
    return sequence_representations[0].tolist()

class LoopData:
    def __init__(self,LP):
        self.__LPD=LP
    def Sequence(self):
        return self.__LPD['Sequence']
    def Length(self):
        return self.__LPD['length']
    def Charge(self):
        return self.__LPD['Charge']
    def Hydro(self):
        return self.__LPD['hydrophobicity']
    def MW(self):
        return self.__LPD['molecular_weight']
    def mz(self):
        return self.__LPD['m/z']
    def Classfication(self):
        return self.__LPD['classification']
    
def BLOSUM62(Loop_Region,idx,DBResult):
    blosum62 = substitution_matrices.load("BLOSUM62")
    Blosum_Scores=[round(pairwise2.align.globalds(Loop_Region,DBS,blosum62,-3,-1)[0].score/len(pairwise2.align.globalds(Loop_Region,DBS,blosum62,-3,-1)[0].seqA),3) for DBS in DBResult]
    PD_Frame=pd.DataFrame({'Loop':DBResult,'Scores':Blosum_Scores})
    PD_Frame.sort_values('Scores', ascending=False)
    Blosum62_Result=PD_Frame.iloc[0:idx]
    return Blosum62_Result,Blosum62_Result['Loop'].values

def NEWSequence(pbar,Blosum_Result,Seq, Clp,BLM):
    Start_index,End_index=Clp.split('-')[0].strip(),Clp.split('-')[1].strip()
    New_Sequence=[Seq[0][:int(Start_index)-1]+LPC+Seq[0][int(End_index):] for LPC in BLM]
    Blosum_Result.insert(0,'Sequence',New_Sequence)
    EMD_Data=[ProteinEmbedding(pbar,[(str(i),New_Sequence[i].strip())]) for i in range(0,len(New_Sequence))]
    return  Blosum_Result,torch.tensor(EMD_Data)

def RunPrediction(Blosum,TT,EMD_Data,Model):
    module=importlib.import_module(f"{Model}.models")
    model1=module.model1
    model1=model1()
    model1.load_state_dict(torch.load(f'{Model}/model1.pth',map_location='cpu',weights_only=True))
    model1.eval()
    model2=module.model2
    model2=model2()
    model2.load_state_dict(torch.load(f'{Model}/model2.pth',map_location='cpu',weights_only=True))
    model2.eval()
    MyModule=module.MyModule
    __main__.MyModule = MyModule
    Pre_Model=torch.load('DAAO/Linear Ensemble.pt',weights_only=False)
    Pre_Model=Pre_Model(model1,model2)
    Blosum['Linear_Regression']=Pre_Model(EMD_Data)
    with open(f'{Model}/Gaussian Ensemble.pkl','rb') as fm:
        gpe=pickle.load(fm)
    gpe_predict=gpe.predict(Pre_Model.RUN(EMD_Data))
    Blosum['Gaussian_Regression']=gpe_predict
    Blosum.to_csv(TT, index=False)

def DBSearch(Seq, Clp,LPD):
    with open(Seq,'r',encoding='utf-8') as ft:
        Sequence=[i.strip() for i in ft.readlines()]
    Start_index,End_index=Clp.split('-')[0].strip(),Clp.split('-')[1].strip()
    Loop_Region=Sequence[0][int(Start_index)-1:int(End_index)-1]
    Loops=LPD.Sequence().values.tolist()
    Length=LPD.Length().values.tolist()
    Classes=LPD.Classfication().values.tolist()
    Loop_class=Loop_Region.strip()[0]+'*'+Loop_Region.strip()[-1]
    Search_Result=[Loops[i]  for i in range(0,LPD.Length().size) if Loop_class==Classes[i] and Length[i] in [len(Loop_Region)-1,len(Loop_Region),len(Loop_Region)+1]]
    return Sequence,Loop_Region,Search_Result

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--Input_Sequence",help="Protein Sequence File")
    parser.add_argument("-c", "--Critical_Loop_Region",help="Critical Loop Region in Input Sequence")
    parser.add_argument("-m", "--Model_type",help="Model Type")
    parser.add_argument("-t", "--Blosum_Loop",type=int,help="Blosum Loops Number")
    parser.add_argument("-o", "--Output_File",help="Predicted File")
    args = parser.parse_args()
    pbar = tqdm(total=100)
    with open('Lp_attr.pkl','rb') as ft:
        Lp_Data=pickle.load(ft)
    Sequence,Lp_Region,DB_Result=DBSearch(args.Input_Sequence,args.Critical_Loop_Region,Lp_Data)
    Blosum_Result,Blosum_Sequence=BLOSUM62(Lp_Region,args.Blosum_Loop,DB_Result)
    pbar.update(int((100-3*args.Blosum_Loop)/3))
    BR,Seq=NEWSequence(pbar,Blosum_Result,Sequence,args.Critical_Loop_Region,Blosum_Sequence)
    pbar.update(int((100-3*args.Blosum_Loop)/3))
    RunPrediction(BR,args.Output_File,Seq,args.Model_type)
    pbar.update(100-3*args.Blosum_Loop-2*int((100-3*args.Blosum_Loop)/3))
    

