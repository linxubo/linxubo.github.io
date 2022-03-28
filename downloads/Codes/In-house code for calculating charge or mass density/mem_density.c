/*
This code mem_density is written using the template "xdrfile-1.1.4.tar.gz" downloaded from www.gromacs.org.
By Xubo Lin (postdoc in Dr. Alex Gorfe's lab, linxbseu@gmail.com)
2017/04/25
Copyright @ Xubo Lin
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xdrfile_xtc.h"

int natoms,step,natom,natom1,read_return;
float time,lambda,f;
matrix box;
rvec *x;
XDRFILE *xtc;
float cDUPC[12],cCLA[1],cSOD[1],cW[1],cDHA[6];            //charge for each atom (unit:e)
float mDUPC[12],mCLA[1],mSOD[1],mW[1],mDHA[6];            //mass for each atom (unit: g/mol)
float tx[150]={0},ty[150]={0},tz[150]={0};                           //1D array to read the coordinates of molecules
float zz[81]={0},ll[81]={0},low[81]={0};                           //Bin width 0.1nm (z range from -4 to 4) for density profiles
int main(int argc,char *argv[])
{
  int mc,nn=0,an,ni,nj,nframes=0;                      //an: atom number; nframes: number of frames used in current analysis
  float t0,t1,bin,ac,am,lmc,sm;             //ac: atom charge; am: atom mass; lmc: COM of lipid bilayer; sm: sum of mass
  FILE *DUPC,*CLA,*SOD,*W,*DHA;
  for(ni=0;ni<81;ni++)
  {
    zz[ni]=-4+0.1*ni;                         //0.1nm: bin width; Set values for output z of density profiles
  }
/*
Read and check the input parameters
*/
  mc=strcmp("mass",argv[3]);                //mc: mass (mc=0) or charge (mc=10)
  if(argc==4 && (mc==0 || mc==10))
  {
    sscanf (argv[1], "%f", &t0);           //change char (e.g. "1234") to int (1234)
    sscanf (argv[2], "%f", &t1);
//  sscanf (argv[3], "%f", &bin);
  }
  else
  {
    printf("%s\n","Please run command to get the mass/charge density: ./mem_density num1 num2 char1");
    printf("%s\n","num1: First frame (ps) to read trajectory; num2: Last frame (ps) to read trajectory; char1: mass or charge");
//  printf("%s\n","Since z output range from -4 to 4 (z=0: center-of-mass of the bilayer), please choose proper bin width!");
  }
/*
Read the charge and mass information of each atom for each molecule
*/
  DUPC=fopen("itp/DUPC.txt","r");
  fscanf(DUPC,"%d %f %f",&an,&ac,&am);
  while(!feof(DUPC))
  {
    cDUPC[nn]=ac,mDUPC[nn]=am;
    nn++;
    fscanf(DUPC,"%d %f %f",&an,&ac,&am);
  }
  fclose(DUPC);
  nn=0;
  W=fopen("itp/W.txt","r");
  fscanf(W,"%d %f %f",&an,&ac,&am);
  while(!feof(W))
  {
    cW[nn]=ac,mW[nn]=am;
    nn++;
    fscanf(W,"%d %f %f",&an,&ac,&am);
  }
  fclose(W);
  nn=0;
  CLA=fopen("itp/CLA.txt","r");
  fscanf(CLA,"%d %f %f",&an,&ac,&am);
  while(!feof(CLA))
  {
    cCLA[nn]=ac,mCLA[nn]=am;
    nn++;
    fscanf(CLA,"%d %f %f",&an,&ac,&am);
  }
  fclose(CLA);
  nn=0;
  SOD=fopen("itp/SOD.txt","r");
  fscanf(SOD,"%d %f %f",&an,&ac,&am);
  while(!feof(SOD))
  {
    cSOD[nn]=ac,mSOD[nn]=am;
    nn++;
    fscanf(SOD,"%d %f %f",&an,&ac,&am);
  }
  fclose(SOD);
  nn=0;
  DHA=fopen("itp/DHA.txt","r");
  fscanf(DHA,"%d %f %f",&an,&ac,&am);
  while(!feof(DHA))
  {
    cDHA[nn]=ac,mDHA[nn]=am;
    nn++;
    fscanf(DHA,"%d %f %f",&an,&ac,&am);
  }
  fclose(DHA);
  nn=0;
/*
Deal with GROMACS xtc trajectory
*/
  nn=0,nframes=0;
  xtc=xdrfile_open("traj.xtc","r");
  read_xtc_natoms("traj.xtc",&natoms);
  x = calloc(natoms, sizeof (x[0]));
  while(1)
  {
    read_return=read_xtc(xtc,natoms,&step,&time,box,x,&f);   //read next frame when using this fuction again. If the current frame is the last frame, the value of read_return is 1
    if(read_return!=0)break;
    if(argc!=4)break;
    if(time<t0)continue;
    if(time>t1)break;
    nframes++;
//    printf("%d\n",nframes);                //test code speed
    lmc=0,sm=0;
/*
Calculate the COM of the lipid bilayer
*/
    for(natom=1;natom<=2760;natom++)             //the bilayer
    {
      lmc=lmc+mDUPC[nn]*x[natom-1][2];              //For density calculation, we only use the z component
      sm=sm+mDUPC[nn];
      nn++;
      if(nn==12)nn=0;
    }
    lmc=lmc/sm;                                  //Center-of-mass of the bilayer (z component)
    sm=0;
/*
Calculate the density profile with COM of the bilayer as origin along z dimension
*/
    for(natom=1;natom<=1332;natom++)             //DUPC
    {
      for(nj=0;nj<81;nj++)
      {
        if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cDUPC[nn];           //charge density
        if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mDUPC[nn];            //mass density
      }
      nn++;
      if(nn==12)nn=0;
    }
    for(natom=1333;natom<=1380;natom++)             //DHA
    {
      for(nj=0;nj<81;nj++)
      {
        if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cDHA[nn];           //charge density
        if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mDHA[nn];            //mass density
      }
      nn++;
      if(nn==6)nn=0;
    }
    for(natom=1381;natom<=2712;natom++)             //DUPC
    {
      for(nj=0;nj<81;nj++)
      {
        if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cDUPC[nn];           //charge density
        if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mDUPC[nn];            //mass density
      }
      nn++;
      if(nn==12)nn=0;
    }
    for(natom=2713;natom<=2760;natom++)             //DHA
    {
      for(nj=0;nj<81;nj++)
      {
        if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cDHA[nn];           //charge density
        if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mDHA[nn];            //mass density
      }
      nn++;
      if(nn==6)nn=0;
    }
    for(natom=2761;natom<=7472;natom++)             //water in the box
    {
      for(nj=0;nj<81;nj++)
      {
        if(abs(x[natom-1][2]-lmc)<4.1)
        {
          if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cW[nn];           //charge density
          if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mW[nn];            //mass density
        }
        if(abs(x[natom-1][2]+box[2][2]-lmc)<4.1)   //In case that [-4,4] crosses the box boundary
        {
          if(mc==10 && (lmc-x[natom-1][2]-box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]-box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cW[nn];           //charge density
          if(mc==0 && (lmc-x[natom-1][2]-box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]-box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mW[nn];            //mass density
        }
        if(abs(x[natom-1][2]-box[2][2]-lmc)<4.1)   //In case that [-4,4] crosses the box boundary
        {
          if(mc==10 && (lmc-x[natom-1][2]+box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]+box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cW[nn];           //charge density
          if(mc==0 && (lmc-x[natom-1][2]+box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]+box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mW[nn];            //mass density
        }
      }
      nn++;
      if(nn==1)nn=0;
    }
    for(natom=7473;natom<=7540;natom++)             //Na+ in the box
    {
      for(nj=0;nj<81;nj++)
      {
        if(abs(x[natom-1][2]-lmc)<4.1)
        {
          if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cSOD[0];           //charge density
          if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mSOD[0];            //mass density
        }
        if(abs(x[natom-1][2]+box[2][2]-lmc)<4.1)  //In case that [-4,4] crosses the box boundary
        {
          if(mc==10 && (lmc-x[natom-1][2]-box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]-box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cSOD[0];           //charge density
          if(mc==0 && (lmc-x[natom-1][2]-box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]-box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mSOD[0];            //mass density
        }
        if(abs(x[natom-1][2]-box[2][2]-lmc)<4.1)  //In case that [-4,4] crosses the box boundary
        {
          if(mc==10 && (lmc-x[natom-1][2]+box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]+box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cSOD[0];           //charge density
          if(mc==0 && (lmc-x[natom-1][2]+box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]+box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mSOD[0];            //mass density
        }
      }
    }
    for(natom=7541;natom<=7592;natom++)             //Cl- in the box
    {
      for(nj=0;nj<81;nj++)
      {
        if(abs(x[natom-1][2]-lmc)<4.1)
        {
          if(mc==10 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cCLA[0];           //charge density
          if(mc==0 && (lmc-x[natom-1][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mCLA[0];            //mass density
        }
        if(abs(x[natom-1][2]+box[2][2]-lmc)<4.1)  //In case that [-4,4] crosses the box boundary
        {
          if(mc==10 && (lmc-x[natom-1][2]-box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]-box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cCLA[0];           //charge density
          if(mc==0 && (lmc-x[natom-1][2]-box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]-box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mCLA[0];            //mass density
        }
        if(abs(x[natom-1][2]-box[2][2]-lmc)<4.1)  //In case that [-4,4] crosses the box boundary
        {
          if(mc==10 && (lmc-x[natom-1][2]+box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]+box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+cCLA[0];           //charge density
          if(mc==0 && (lmc-x[natom-1][2]+box[2][2])<=zz[nj]+0.025 && (lmc-x[natom-1][2]+box[2][2])>=zz[nj]-0.025)ll[nj]=ll[nj]+mCLA[0];            //mass density
        }
      }
    }
    for(ni=0;ni<81;ni++)
    {
      low[ni]=low[ni]+ll[ni]/(0.1*box[0][0]*box[1][1]);        //0.1: bin width
      ll[ni]=0;
    }
  }
  xdrfile_close(xtc);
/*
Output the density profile with the format: z coordiantes (-4 to 4), charge/mass density of the lipid bilayer
*/
  for(ni=0;ni<81;ni++)
  {
    if(mc==10)printf("%f %f\n",zz[ni],low[ni]/nframes);        //charge density profiles, unit: e/nm3
    if(mc==0)printf("%f %f\n",zz[ni],10*low[ni]/(6.02214*nframes));        //mass density profiles, unit: kg/m3
  }
}
