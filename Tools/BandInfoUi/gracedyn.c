#include <grace_np.h>
#include <stdio.h>

int DGraceOpen(int bs)
{
  return(GraceOpen(bs));
}


int DGraceOpen2(int bs)
{
//  return(GraceOpenVA("./xx",bs,"-noask", NULL));
}


int DGraceCommand(const char* cmd){
//  printf("%s\n",cmd);
//  fflush(stdout);
  return(GraceCommand(cmd));  
}
