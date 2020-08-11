/*#include "ProdTutorial/ProducerTest/plugins/DetImgProducer.h"
#include "ProdTutorial/ProducerTest/plugins/EGProducer.h"
#include "ProdTutorial/ProducerTest/plugins/QGProducer.h"*/
#include "ProdTutorial/ProducerTest/plugins/EGProducer.h"
#include "ProdTutorial/ProducerTest/plugins/QGProducer.h"
//using namespace std;

std::vector<std::vector<float>> croppingFrames(std::vector<float>& vdetector_image, int ieta_seed, int iphi_seed, int detImg_height, int detImg_width, int frame_height, int frame_width){
  std::vector<std::vector<float>> vframe = std::vector<std::vector<float>> (frame_height,std::vector<float> (frame_width, 0.0));
  int start_x=0;
  int end_x=0;
  int start_y=0;
  //int end_y=0;
  int buff_x=0;
  int buff_y=0;
  int half_frame_width;
  int half_frame_height;
  if (frame_width%2==0){half_frame_width=frame_width/2-1;} else {half_frame_width=frame_width/2;}
  if(frame_height%2==0){half_frame_height=frame_height/2-1;} else {half_frame_height=frame_height/2;}
  if (iphi_seed<half_frame_width-1){
   start_y=0;
   buff_y=half_frame_width-iphi_seed;
   buff_y=detImg_width-buff_y;
  }
  else {
   start_y=iphi_seed-half_frame_width;
   buff_y=0;
  }
  if (iphi_seed>detImg_width-half_frame_width){
   //end_y=detImg_width-1;
  }
  else {
   //end_y=iphi_seed+frame_width/2;
  }
  if (ieta_seed<half_frame_height){
   start_x=0;
   buff_x=half_frame_height-ieta_seed;
  }
  else {
   start_x=ieta_seed-half_frame_height;
    buff_x=0;
  }
  if (ieta_seed>detImg_height-half_frame_height){
   end_x=detImg_height-1;
  }
  else {
   end_x=ieta_seed+half_frame_height;
  }
  /*std::string filename = "frame_" + std::to_string(iP+1) + "_" + std::to_string(nPassed+1) + ".csv";
  std::ofstream frame_file(filename);*/
  for (int x_idx = start_x; x_idx<=end_x;x_idx++){
   for (int y_idx = 0/*start_y*/; y_idx<frame_width/*=end_y*/;y_idx++){
    vframe[x_idx-start_x+buff_x][y_idx/*y_idx-start_y+buff_y*/]=vdetector_image[x_idx*detImg_width+(y_idx+buff_y+start_y)%detImg_width];
    //vEB_flat_frame[(x_idx-start_x+buff_x)*frame_width+y_idx/*-start_y+buff_y*/]=vdetector_image[x_idx*detImg_width+(y_idx+start_y+buff_y)%detImg_width];
    //std::cout<<"("<<x_idx-start_x+buff_x<<","<<y_idx<<"): "<<vframe[x_idx-start_x+buff_x][y_idx/*y_idx-start_y+buff_y*/]<<" "<<vdetector_image[x_idx*detImg_width+(y_idx+start_y+buff_y)%detImg_width];
   }
   //std::cout<<std::endl;
  }
  /*for (int x_idx=0;x_idx<frame_height;x_idx++){
   for (int y_idx=0;y_idx<frame_width;y_idx++){ 
    std::cout<<"("<<x_idx<<","<<y_idx<<"): "<<vframe[x_idx][y_idx]<<" ";
    frame_file<<vframe[x_idx][y_idx];
    if (y_idx<(frame_width-1)){frame_file<<",";}
   }
    frame_file<<"\n";
  }*/
  //vEB_photon_frames.push_back(vEB_flat_frame);
  std::cout<<" >> Size of frame is:"<<"("<<vframe.size()<<", "<<vframe[0].size()<<")"<<endl;
  //std::cout<<" >> E_max at ("<<ieta_seed<<", "<<iphi_seed<<")is: "<<vdetector_image[ieta_seed*detImg_width+iphi_seed]<<std::endl;
  std::cout<<" >> Detector Image value at ("<<ieta_seed<<", "<<iphi_seed<<")is: "<<vframe[half_frame_height][half_frame_width]<<std::endl;
  //std::cout<<ieta_seed<<" "<<iphi_seed<<" "<<buff_y<<" "<<start_y<<std::endl;
  return vframe;
}

