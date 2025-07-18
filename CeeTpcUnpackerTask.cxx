//-----------------------------------------------------------
// File and Version Information:
// $Id$

// Description:
//      Implementation of class TpcClusterFinderTask
//
// Environment:
//      Software developed for the MPD at NICA.
//
// Author List:
//      Sergey Merts
//      Artem Basalaev
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "CeeDAQ_BasicUnpacker.h"
#include "CeeTpcUnpackerTask.h"
#include "CeeTpcDigi.h"
//#include "TpcCluster.h"
#include "FairRootManager.h"
#include "CeeMCTrack.h"
#include "TClonesArray.h"
#include <TSystem.h>
#include <iostream>
#include "TaskHelpers.h"
#include <TGraph.h>
#include "TMath.h"
#include "TVector3.h"
#include "McIdCollection.h"
#include "McId.h"
#include "sys/timeb.h"
#include <sys/time.h>
#include <fstream>
#include <ostream>
#include "TStopwatch.h"
#include "TString.h"
#include <ios>
using namespace std;
using namespace TMath;
using namespace CeeTpcUnpackerNameSpace;

ClassImp(CeeTpcUnpackerTask);


CeeTpcUnpackerTask::CeeTpcUnpackerTask() :
fPersistence(kTRUE),
fPrintDebugInfo(kFALSE),
fDigisArray(NULL),
isHistogramsInitialized(kFALSE),
fMakeQA(kFALSE),
is(NULL)
{
    outputBranchName = "CeeTpcDigi";
    inputBranchName  = "CeeT0Digi";
    map_infile = gSystem->Getenv("VMCWORKDIR");
    // map_infile += "/geometry/tpc_map_1760.txt";
    map_infile += "/tpc/unpacker/ROUNDPAD15312_MapFile_Sorted_by_Electronics_v1.2.root";
}

CeeTpcUnpackerTask::~CeeTpcUnpackerTask() {
    // if (isHistogramsInitialized) delete fHisto;
    delete fDigisArray;
    // delete fGas;
}

InitStatus CeeTpcUnpackerTask::Init() 
{
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman) 
    {
        cout << "\n-E- [CeeTpcUnpackerTask::Init]: RootManager not instantiated!" << endl;
        return kFATAL;
    }

    fT0Array = (TClonesArray*) ioman->GetObject(inputBranchName);
    fDigisArray = new TClonesArray(outputBranchName);
    ioman->Register(outputBranchName, "CeeTpc", fDigisArray, fPersistence);
    fEventId = 0;
    ReadMapFile(map_infile);
//=========================================================================================================================================================================================

///// Event Header >>>>>
//=========================================================================================================================================================================================

    cout << "-I- CeeTpcUnpackerTask: Initialization successful." << endl;
    return kSUCCESS;
}
void CeeTpcUnpackerTask::Exec(Option_t* opt) {
    cout << "CeeTpcUnpacker::Exec started! Event #" << fEventId << endl;
    if (!fDigisArray) Fatal("CeeTpcUnpacker::Exec)", "No FoundClustersArray");
    fDigisArray->Delete();   
//=========================================================================================================================================================================================
////// Get T0-leading time from CeeT0Digi    
    Double_t Sum_Tleading=0;
    for (register int i = 0; i < fT0Array->GetEntriesFast(); ++i) 
    {
       CeeT0Digi* fT0Digit = (CeeT0Digi*) fT0Array->At(i);
       Sum_Tleading += fT0Digit->GetTime();
    }
    Double_t Mean_Tleading = Sum_Tleading/(fT0Array->GetEntriesFast()); // [ns]
////// TPC data uppacking <<<<<<<<<<<<<<<<
    {
        Int_t event_index_of_all = basic_unpacker->GetEventIndexofAll();
        int file_index = basic_unpacker->GetFileIndex();
        is = basic_unpacker->GetStream();
        if (!is->is_open()) Fatal("CeeTpcUnpacker::Exec)", "No FoundClustersArray");
        CeeDAQ_EventHead* event_header = basic_unpacker->GetEventHeader();
        Int_t event_index = basic_unpacker->GetEventIndex();
        fEventId++;
        uint64_t now_triggerid = event_header->trigerid;
        int child_package_head_length = 8; // TPC数据净荷::子包包头长度（8 byte）
        int child_package_subtriggerid_length = 2; // TPC数据净荷::子包子触发号长度(2 byte)
        AssembleCollectionHead this_TPC_Assemble;
        // {
        //     auto TPC_Assemble_1 = event_header->assemble_list.find(assemble_id1);
        //     auto TPC_Assemble_2 = event_header->assemble_list.find(assemble_id2);
        //     if(TPC_Assemble_1 == event_header->assemble_list.end() && TPC_Assemble_2==event_header->assemble_list.end()){
        //         Info("CeeTpc_Unpacker::Exec","No 0x30a 0x30b trigger%lu\n",now_triggerid);
        //         if(fPrintDebugInfo) Info("CeeTpc_Unpacker::Exec","No TPC Data in Event %d of File %d\n",event_index,file_index);
        //         return;
        //     }
        //     if(TPC_Assemble_1!=event_header->assemble_list.end())
        //     this_TPC_Assemble.Add(TPC_Assemble_1->second);
        //     if(TPC_Assemble_2!=event_header->assemble_list.end())
        //     this_TPC_Assemble.Add(TPC_Assemble_2->second);
        // }
        for(const auto& iassid : assemble_id_list){
            auto i_assemble = event_header->assemble_list.find(iassid);
            if (i_assemble != event_header->assemble_list.end()){
                this_TPC_Assemble.Add(i_assemble->second);
            }
            else {
                if(fPrintDebugInfo)
                printf("TPC: No %#lX in trigger%lu\n",iassid,now_triggerid);
            }
        }
        auto FEE_num = this_TPC_Assemble.FEE_num; // 此触发下tpc的有效fee数量
        for(int FEE_index=0;FEE_index<FEE_num;FEE_index++){ // 遍历FEE
            int BDM_id = this_TPC_Assemble.FEE_id[FEE_index];
            uint64_t FEE_address = this_TPC_Assemble.FEE_address[FEE_index];
            uint64_t FEE_data_length = this_TPC_Assemble.FEE_data_length[FEE_index];
            uint64_t valid_ch_num=0;
            uint64_t fee_data_length=0;
            unsigned char child_package_head[child_package_head_length]; // TPC数据净荷::子包包头（8 byte）
            unsigned char child_package_subtriggerid[child_package_subtriggerid_length]; // TPC数据净荷::子包子触发号（2 byte）
            ULong64_t child_package_6ByteZero=0;
            uint64_t RealData_address = FEE_address+2; // 跳过FEE板号（2 byte）和FEE数据段长度（3 byte）到达数据净荷位置
            is->seekg(RealData_address, is->beg);
            is->read((char*)&fee_data_length,3); // FEE数据段长度（3 byte）
            is->read((char*)&valid_ch_num,2); // TPC数据净荷::通道数（2 byte）
            if(fPrintDebugInfo) cout << Form("In Fee_%d, valid_ch_num: %lu, fee_data_length: %lu", BDM_id, valid_ch_num, fee_data_length) << endl;
            if(valid_ch_num==0) continue;
            uint64_t child_package_position = RealData_address+2;
            auto FEE_trigger_time = this_TPC_Assemble.FEE_trigertime.at(FEE_index);
            // T0_leading_time to TPC_FEE trigger time
            Double_t T0leading_TPCFEE = FEE_trigger_time*25 - Mean_Tleading; // [ns]
            if(fPrintDebugInfo) cout << Form(" TPCFEE-T0: %e ,FEE_trigger_time: %e ,Mean_T0 leading: %e",T0leading_TPCFEE,FEE_trigger_time,Mean_Tleading) << endl;
            // T0 to SAMPA
            Double_t T0_SAMPA = T0leading_TPCFEE + 93.; // [ns]
            if(!fUseT0) T0_SAMPA = 0;
            // TPC sample latency(ns)
            Double_t TPC_sample_latency = -1*192.*100.; // [ns]
            for(int ch_index=0;ch_index<valid_ch_num;ch_index++){ // 遍历子包
                // is->seekg(child_package_position, is->beg);
                // cout << is->tellg() << endl;
                is->read((char*)(&child_package_head), child_package_head_length); // 读取子包包头（8 byte）
                is->read((char*)(&child_package_subtriggerid), child_package_subtriggerid_length); // 读取子包子触发号（2 byte）
                is->read((char*)&child_package_6ByteZero,6);// 保留字节 6 byte
                if(child_package_head[7] != 0x40 || child_package_6ByteZero!=0) {
                    if(fPrintDebugInfo) cout << Form("Not child_package_head: %x", child_package_head[7]) << endl;
                    break;
                }
                // printf("Is child_package_head: %x\n",child_package_head[7]);
                DAQ_TPC_Channel temp_channel;
                UInt_t subtriggerid = 0;
                memcpy(&subtriggerid, child_package_subtriggerid, 2);
                temp_channel.SubTrigger_ID = subtriggerid;
                temp_channel.Trigger_ID = now_triggerid;
                temp_channel.Fee_ID = BDM_id;
                temp_channel.Chip_ID = (child_package_head[2] & 0xf0) >> 4;
                temp_channel.Channel_ID = child_package_head[3] & 0x1f;
                temp_channel.BX_counter = ((child_package_head[3]&0x20)>>5) + ((child_package_head[4]&0xff)<<1) + ((child_package_head[5]&0xff)<<9) + ((child_package_head[6]&0x07)<<17);
                Double_t T0_TPCWindow = (T0_SAMPA + temp_channel.BX_counter*25 + TPC_sample_latency)*1e-2; // [100ns] same as Sample_timebin
                int sci_data_length = ((child_package_head[1]&0xfc) >> 2) + ((child_package_head[2]&0x0f) << 6); // 科学数据长度，以10 bit为单位，6个10 bit对应8 byte
                int quotient = sci_data_length / 6;
                int remainder = sci_data_length % 6;
                int real_sci_data_lenght = 0; // 以byte为单位的科学数据长度，8的整数倍
                real_sci_data_lenght = (quotient + ((remainder > 0)?1:0)) * 8;
                unsigned char sci_data[real_sci_data_lenght];
                is->read((char*)(&sci_data), real_sci_data_lenght);
                for(int sci_data_index = 0; sci_data_index < real_sci_data_lenght; sci_data_index += 8){
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                    temp_channel.sci_data.push_back((sci_data[sci_data_index+0] & 0xff) + ((sci_data[sci_data_index+1] & 0x03)<<8));
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                    temp_channel.sci_data.push_back(((sci_data[sci_data_index+1] & 0xfc)>>2) + ((sci_data[sci_data_index+2] & 0x0f)<<6));
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                    temp_channel.sci_data.push_back(((sci_data[sci_data_index+2] & 0xf0)>>4) + ((sci_data[sci_data_index+3] & 0x3f)<<4));
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                    temp_channel.sci_data.push_back((sci_data[sci_data_index+4] & 0xff) + ((sci_data[sci_data_index+5] & 0x03)<<8));
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                    temp_channel.sci_data.push_back(((sci_data[sci_data_index+5] & 0xfc)>>2) + ((sci_data[sci_data_index+6] & 0x0f)<<6));
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                    temp_channel.sci_data.push_back(((sci_data[sci_data_index+6] & 0xf0)>>4) + ((sci_data[sci_data_index+7] & 0x3f)<<4));
                    if(temp_channel.sci_data.size() == sci_data_length) break;
                }
                bool is_badchild = 0;
                {
                    if(temp_channel.sci_data.size()<3) is_badchild=true; // 科学数据总长<3时，显然是坏数据
                    else{
                        int i_scidata=0;
                        unsigned int cluster_n=0;
                        unsigned int cluster_t=0;
                        while (i_scidata<temp_channel.sci_data.size())
                        {
                            cluster_n = temp_channel.sci_data.at(i_scidata);
                            if(cluster_n==0){
                                is_badchild=true; // cluster长度为0,显然是坏数据
                                break;
                            }
                            else{
                                i_scidata+=(1+cluster_n);
                                if(i_scidata>=temp_channel.sci_data.size()){
                                    is_badchild=true; // 当前cluster最后一个元素位置溢出,显然是坏数据
                                    break;
                                }
                                else i_scidata++;
                            }
                        }
                    }
                }
                if(is_badchild) continue;
                if(fPrintDebugInfo) std::cout << Form("ch:%u-%u-%u\tchildhead:",temp_channel.Fee_ID,temp_channel.Chip_ID,temp_channel.Channel_ID);
                for(int i_child=0;i_child<child_package_head_length;i_child++){
                    if(fPrintDebugInfo) std::cout << Form("%02X",child_package_head[i_child]);
                }
                if(fPrintDebugInfo) std::cout << std::endl << "sci_data: ";
                for(auto i_scidata=temp_channel.sci_data.begin();i_scidata!=temp_channel.sci_data.end();++i_scidata){
                    if(fPrintDebugInfo) std::cout << Form("%u ",*i_scidata);
                }
                if(fPrintDebugInfo) std::cout << std::endl;
                // Get CeeTpcDigi from temp_channal
                // ALineInMapFile mapline;
                // if(temp_channel.GetAMapLine(map_infile,mapline))
                auto mapline = GetMapLine(temp_channel.Fee_ID*128 + temp_channel.Chip_ID*32 + temp_channel.Channel_ID);
                if(mapline.TH2POLY_BIN_ID > 0) {
                    for(int sci_data_index = 0; sci_data_index < temp_channel.sci_data.size(); )
                    { // loop sci_data
                        auto cluster_num = temp_channel.sci_data[sci_data_index];
                        auto cluster_t = temp_channel.sci_data[sci_data_index+1];
                        auto now_sample_index = sci_data_index + 2;
                        auto next_index = sci_data_index + 2 + cluster_num;
                        for(sci_data_index += 2; sci_data_index < next_index; ++sci_data_index){ // loop in a cluster
                            Int_t Sample_timebin = cluster_t + sci_data_index - now_sample_index;
                            Int_t Sample_adc = temp_channel.sci_data[sci_data_index];
                            CeeTpcDigi* outDigi;
                            {
                                // outDigi = new((*fDigisArray)[fDigisArray->GetEntriesFast()]) CeeTpcDigi(mapline.BinID-1, mapline.IDInRow-1, 0, mapline.RowNum-1,Sample_timebin+T0_TPCWindow, Sample_adc, 0);
                                outDigi = new((*fDigisArray)[fDigisArray->GetEntriesFast()]) CeeTpcDigi(mapline.TH2POLY_BIN_ID-1, mapline.col_id, mapline.sec_id, mapline.row_id, Sample_timebin+T0_TPCWindow, Sample_adc, 0);
                                if(fPrintDebugInfo) std::cout<<"mapline.TH2POLY_BIN_ID = "<<mapline.TH2POLY_BIN_ID-1<< "    Sample_timebin = "<<Sample_timebin<<"    Sample_adc = "<<Sample_adc<< std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
////// TPC data uppacking end <<<<<<<<<<<<<<<<
//=========================================================================================================================================================================================

    // Reset output Array
    if (fPrintDebugInfo) cout << "Number of Found digi = " << fDigisArray->GetEntriesFast() << endl;
}


//===========================================================================================================================================================================================================
void CeeTpcUnpackerTask::Finish() {
cout<<"unpacker successed"<<endl;
    
//mybeg
/*
    TH2F *_thist = new TH2F("thist","thist_title",300,-150,150,300,-150,150); 
    for(int i=0; i<fDigisArray->GetEntriesFast(); i++) {
        CeeTpcHit *fHit = (CeeTpcHit*) fDigisArray->At(i); 
        _thist->Fill(fHit->GetX(), fHit->GetY(), fHit->GetQ());
        tfile2 <<i<<") "<<"XYZ="<<fHit->GetLocalX()<<", "<<fHit->GetLocalY()<<", "<<fHit->GetLocalZ()<<";  Layer= "<<fHit->GetLayer()<< ";  Pad="<<fHit->GetPad()<<";  Bin"<<fHit->GetBin()<<"\n";
    }
    
    _thist->Write("THIST", kOverwrite);

*/ 
//myend    
}

void CeeTpcUnpackerTask::ReadMapFile(const TString& path_mapfile){
    // std::unique_ptr<TFile> map_file( TFile::Open(path_mapfile.Data(), "READ") );
    // auto maptree = map_file->Get<TTree>("maptree");
    auto maptree = new TChain("maptree");
    maptree->Add(path_mapfile.Data());
    PadInfo mappad;
    maptree->SetBranchAddress("sector",   &mappad.sec_id);
    maptree->SetBranchAddress("row",      &mappad.row_id);
    maptree->SetBranchAddress("column",   &mappad.col_id);
    maptree->SetBranchAddress("th2polybin",   &mappad.TH2POLY_BIN_ID);
    maptree->SetBranchAddress("fee",      &mappad.FEE_ID);
    maptree->SetBranchAddress("sampa",   &mappad.SAMPA_ID);
    maptree->SetBranchAddress("channel",      &mappad.CHANNEL_ID);
    // maptree->SetBranchAddress("pm_x",     &mappad.pm.x);
    // maptree->SetBranchAddress("pm_z",     &mappad.pm.z);
    map_sorted.clear();
    for(int ipad=0;ipad<maptree->GetEntries();ipad++){
        maptree->GetEntry(ipad);
        auto tmp_pad = mappad;
        map_sorted.push_back(tmp_pad);
    }
}
