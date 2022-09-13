from PyQt5 import QtWidgets, uic

from PyQt5.QtWidgets import QFileDialog

import sys, os

import ctypes



class Ui(QtWidgets.QMainWindow):

    def close(self):
      print ("closing")
      self.lib.GraceClose()
      super(Ui, self).close()

    def __init__(self):


        self.cwd = os.getcwd();


        self.lib = ctypes.CDLL("./gracedyn.so")
        self.lib.GraceOpen(2048);

        super(Ui, self).__init__()
        uic.loadUi('BandInfo.ui', self)

        self.setGeometry(900,0,self.width(),self.height())

        self.show()

        self.all_S =0
        self.all_P =0
        self.all_D =0
        self.all_F =0

        self.do_S = 0
        self.do_P = 0
        self.do_D = 0
        self.do_F = 0


        self.do_S0  = 0

        self.do_P_1 = 0
        self.do_P0  = 0
        self.do_P1  = 0

        self.do_D_2  = 0
        self.do_D_1  = 0
        self.do_D0   = 0
        self.do_D1   = 0
        self.do_D2   = 0

        self.do_F_3  = 0
        self.do_F_2  = 0
        self.do_F_1  = 0
        self.do_F0   = 0
        self.do_F1   = 0
        self.do_F2   = 0
        self.do_F3   = 0

        self.do_SpinDown=1
        self.do_SpinUp=1

        self.checkBoxSpec1.hide()
        self.checkBoxSpec2.hide()
        self.checkBoxSpec3.hide()
        self.checkBoxSpec4.hide()
        self.checkBoxSpec5.hide()
        self.checkBoxSpec6.hide()
        self.checkBoxSpec7.hide()
        self.checkBoxSpec8.hide()


    #    self.InitList()


    def genLineEditSpecies(self):
      set_str=""
      if (self.checkBoxSpec1.isChecked()) :
        set_str = set_str + "1,"
      if (self.checkBoxSpec2.isChecked()) :
        set_str = set_str + "2,"
      if (self.checkBoxSpec3.isChecked()) :
        set_str = set_str + "3,"
      if (self.checkBoxSpec4.isChecked()) :
        set_str = set_str + "4,"
      if (self.checkBoxSpec5.isChecked()) :
        set_str = set_str + "5,"
      if (self.checkBoxSpec6.isChecked()) :
        set_str = set_str + "6,"
      if (self.checkBoxSpec7.isChecked()) :
        set_str = set_str + "7,"
      if (self.checkBoxSpec8.isChecked()) :
        set_str = set_str + "8,"


      set_str = set_str[:-1]
      self.lineEditSpecies.setText(set_str)


    def InitList(self):

        if (self.checkBoxSpin.isChecked()) :
          self.pushButtonSpinUp.setDown(True)
          self.pushButtonSpinDown.setDown(True)
        else:
          self.pushButtonSpinUp.setDown(False)
          self.pushButtonSpinDown.setDown(False)


        fquery = open(self.dirname+"/"+"band_info_query.dat","w")

        if (self.checkBoxSpin.isChecked()) :
          fquery.write("2\n")
          fquery.write("n\n")
        else:
          fquery.write("1\n")
          fquery.write("n\n")

        fquery.close()

        os.system("cd " + self.dirname + "; " + self.cwd+"/BandInfo.exe < band_info_query.dat ")

        fin = open(self.dirname+"/"+"band_info.txt","r")

        self.neig =  int(fin.readline());
        self.nbaslcao = int(fin.readline());

        print ("neig", self.neig, "nbaslcao", self.nbaslcao)

        self.info_itype = [None]*self.nbaslcao
        self.info_iatom = [None]*self.nbaslcao
        self.info_n     = [None]*self.nbaslcao
        self.info_l     = [None]*self.nbaslcao
        self.info_m     = [None]*self.nbaslcao

        for iorb in range(0,self.nbaslcao):
          info = fin.readline().split();
          self.info_itype[iorb]=int(info[1])
          self.info_iatom[iorb]=int(info[2])
          self.info_n[iorb]=int(info[3])
          self.info_l[iorb]=int(info[4])
          self.info_m[iorb]=int(info[5])


        self.nspecies = int(fin.readline())

        self.nameat =[None]*(self.nspecies+1)

        self.name_keys = {}

        for ispecies in range(0,self.nspecies):
          self.nameat[ispecies] = fin.readline().rstrip()
          self.name_keys[self.nameat[ispecies]] = ispecies


        fin.close()

#        print ("name_keys=", name_keys)
#        exit()


        for ispecies in range(0,self.nspecies):
          if (ispecies == 0) :
            self.checkBoxSpec1.show()
            self.checkBoxSpec1.setText(self.nameat[ispecies])
            self.checkBoxSpec1.setChecked(True)
          if (ispecies == 1) :
            self.checkBoxSpec2.show()
            self.checkBoxSpec2.setText(self.nameat[ispecies])
            self.checkBoxSpec2.setChecked(True)
          if (ispecies == 2) :
            self.checkBoxSpec3.show()
            self.checkBoxSpec3.setText(self.nameat[ispecies])
            self.checkBoxSpec3.setChecked(True)
          if (ispecies == 3) :
            self.checkBoxSpec4.show()
            self.checkBoxSpec4.setText(self.nameat[ispecies])
            self.checkBoxSpec4.setChecked(True)
          if (ispecies == 4) :
            self.checkBoxSpec5.show()
            self.checkBoxSpec5.setText(self.nameat[ispecies])
            self.checkBoxSpec5.setChecked(True)
          if (ispecies == 5) :
            self.checkBoxSpec6.show()
            self.checkBoxSpec6.setText(self.nameat[ispecies])
            self.checkBoxSpec6.setChecked(True)
          if (ispecies == 6) :
            self.checkBoxSpec7.show()
            self.checkBoxSpec7.setText(self.nameat[ispecies])
            self.checkBoxSpec7.setChecked(True)
          if (ispecies == 7) :
            self.checkBoxSpec8.show()
            self.checkBoxSpec8.setText(self.nameat[ispecies])
            self.checkBoxSpec8.setChecked(True)



        print("self.nspecies is ", self.nspecies)

        print("self.nameat is ", self.nameat)

        natom=[None]*(self.nspecies+1)


        for ispecies in range(1,self.nspecies+1):
          natom[ispecies]=0
          print("ispecies natom", ispecies,natom)
          if (ispecies ==1):
            self.lineEditSpecies.setText( str(ispecies))
          if (ispecies >1):
            self.lineEditSpecies.setText( self.lineEditSpecies.text() +"," + str(ispecies))

        for iorb in range(0,self.nbaslcao):
          ispecies = self.info_itype[iorb]
          if (self.info_iatom[iorb] > natom[ispecies]) :
            natom[ispecies]=natom[ispecies] +1


        atom_list_str = "";

        for ispecies in range(1, self.nspecies+1):
          for iatom in range(1, natom[ispecies]+1):
#            print ("iatom=", iatom)
            atom_list_str = atom_list_str + str(iatom) + " "
            if (iatom < natom[ispecies]):
              atom_list_str = atom_list_str + ","
            if (iatom == natom[ispecies]):
              atom_list_str = atom_list_str + "\n"


#        print ("atom_list_str", atom_list_str)

        self.textEditAtoms.setPlainText(atom_list_str)

#        print("natom list", natom)

        self.lineEditSpecies.setEnabled(False)
        self.textEditAtoms.setEnabled(False)


    def CheckAtoms(self,icheck):
      if (icheck ==0 ) :
        self.textEditAtoms.setEnabled(True)
      else:
        self.textEditAtoms.setEnabled(False)

    def CheckSpecies(self,icheck):
      if (icheck ==0 ) :
        self.lineEditSpecies.setEnabled(True)
      else:
        self.lineEditSpecies.setEnabled(False)

    def Load(self):

        self.dirname = QFileDialog.getExistingDirectory()

        self.InitList()

        self.pushButtonRun.setEnabled(True)

    def Run(self):

        iselect=[None]*self.nbaslcao
        for iorb in range(0,self.nbaslcao):
          iselect[iorb]=0


        SpeciesList  = self.lineEditSpecies.text().split(",");

        print ("SpecList", SpeciesList, len(SpeciesList))

        if(SpeciesList[0]==""):
            return

#        print ("Species List is", SpeciesList)



        AllAtomList  = self.textEditAtoms.toPlainText().split("\n");
#        print ("Atom List is", AllAtomList)


        for iorb in range(0,self.nbaslcao):
          for species in SpeciesList:
            ispecies = int(species)
            AtomList = AllAtomList[ispecies-1].split(",")
#            print ("iatom_list", AtomList )
            for atom in AtomList:
              iatom = int(atom)
#              print("iatom", iatom)
              if(self.info_itype[iorb] == ispecies):
                if(self.info_iatom[iorb] == iatom):
                  if(self.info_l[iorb] == 0) :
                    if(self.do_S0 == 1):
                      iselect[iorb] = 1
                  if(self.info_l[iorb] == 1) :
                    if(self.info_m[iorb] == -1 and self.do_P_1 ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] ==  0 and self.do_P0  ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] ==  1 and self.do_P1  ==1):
                      iselect[iorb] = 1
                  if(self.info_l[iorb] == 2) :
                    if(self.info_m[iorb] == -2 and self.do_D_2 ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] == -1 and self.do_D_1 ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] == 0  and self.do_D0  ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] == 1  and self.do_D1  ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] == 2  and self.do_D2  ==1):
                      iselect[iorb] = 1
                  if(self.info_l[iorb] == 3) :
                    if(self.info_m[iorb] == -3 and self.do_F_3 ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] == -2 and self.do_F_2 ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] == -1 and self.do_F_1 ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] ==  0 and self.do_F0  ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] ==  1 and self.do_F1  ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] ==  2 and self.do_F2  ==1):
                      iselect[iorb] = 1
                    if(self.info_m[iorb] ==  3 and self.do_F3  ==1):
                      iselect[iorb] = 1


        if (self.checkBoxSpin.isChecked()) :
          for iorb in range(0,int(self.nbaslcao/2)):
            print ("irob", iorb, self.do_SpinUp, self.do_SpinDown)
            if (self.do_SpinUp==0) :
              iselect[2*iorb]=0
            if (self.do_SpinDown==0) :
              iselect[2*iorb+1]=0

        for iorb in range(0,self.nbaslcao):
          print (iorb, "itype", self.info_itype[iorb], "iatom", self.info_iatom[iorb], "n", self.info_n[iorb], "l", self.info_l[iorb], "m", self.info_m[iorb], "iselect", iselect[iorb] )

        fout = open(self.dirname+"/"+"band_info_run.dat", "w")

        if (self.checkBoxSpin.isChecked()) :
          fout.write("2\n")
        else:
          fout.write("1\n")

        fout.write("y\n")
        fout.write("y\n")

        for iorb in range(0,self.nbaslcao):
          fout.write(str(iselect[iorb])+"\n")

        fout.close()


        os.system("cp " + self.cwd + "/Default.agr " + self.dirname +  " ; cd " + self.dirname + "; pwd ; " + self.cwd+"/BandInfo.exe < band_info_run.dat > b_log");
        if (self.checkBoxSpin.isChecked()) :
          strload = "load \""+ self.dirname+"/"+"band_fld_so_all.agr\""
          self.lib.GraceCommand(str.encode(strload));
        else:
          strload = "load \"" + self.dirname+"/"+"band_fld_all.agr\""
          self.lib.GraceCommand(str.encode(strload));

        self.lib.GraceCommand(b"redraw");

        return


#-----
        # f = open("rep.dat", "w")
        # f.write("   1    Which file do you want to read ?\n")
        # f.write("   y    do you want to proceed? ?\n")

        # f.write("y   do you want to see all atoms?\n")

        # if(self.do_S == 1) :
          # f.write("   0    l quantum number (-1 for all)\n")
          # f.write(" -10    m quantum number (-10 for all)\n")
        # elif(self.do_P == 1) :
          # f.write("   1    l quantum number (-1 for all)\n")
          # f.write(" -10    m quantum number (-10 for all)\n")
        # elif(self.do_D == 1) :
          # f.write("   2    l quantum number (-1 for all)\n")
          # f.write(" -10    m quantum number (-10 for all)\n")


        # f.close()
        # self.CountM()
#        return


    def CountM(self):
        ans = self.do_P_1 and self.do_P0 and self.do_P1
        if (ans) :
          self.all_P = True
        else:
          self.all_P = False

        ans = self.do_D_2 and self.do_D_1 and self.do_D0 and self.do_D1 and self.do_D2
        if (ans) :
          self.all_D = True
        else:
          self.all_D = False

        if (self.do_S0) :
          self.all_S = True
        else:
          self.all_S = False


        print("ans ")
        print(self.all_S)
        print("ans ")
        print(self.all_P)
        print("ans ")
        print(self.all_D)
        print("ans ")
        print(self.all_F)


    def DisableAllLandM(self):

        self.do_S = 0
        self.do_P = 0
        self.do_D = 0
        self.do_F = 0


        self.do_S0  = 0

        self.do_P_1 = 0
        self.do_P0  = 0
        self.do_P1  = 0

        self.do_D_2  = 0
        self.do_D_1  = 0
        self.do_D0   = 0
        self.do_D1   = 0
        self.do_D2   = 0

        self.do_F_3  = 0
        self.do_F_2  = 0
        self.do_F_1  = 0
        self.do_F0   = 0
        self.do_F1   = 0
        self.do_F2   = 0
        self.do_F3   = 0


        self.pushButtonS.setDown(False)
        self.pushButtonS0.setDown(False)

        self.pushButtonP.setDown(False)
        self.pushButtonP_1.setDown(False)
        self.pushButtonP0.setDown(False)
        self.pushButtonP1.setDown(False)

        self.pushButtonD.setDown(False)
        self.pushButtonD_2.setDown(False)
        self.pushButtonD_1.setDown(False)
        self.pushButtonD0.setDown(False)
        self.pushButtonD1.setDown(False)
        self.pushButtonD2.setDown(False)

        self.pushButtonF.setDown(False)
        self.pushButtonF_3.setDown(False)
        self.pushButtonF_2.setDown(False)
        self.pushButtonF_1.setDown(False)
        self.pushButtonF0.setDown(False)
        self.pushButtonF1.setDown(False)
        self.pushButtonF2.setDown(False)
        self.pushButtonF3.setDown(False)

    def SpinDown_clicked(self):
        if ( self.do_SpinDown == 0):
          self.pushButtonSpinDown.setDown(True)
          self.do_SpinDown = 1
        else:
          self.pushButtonSpinDown.setDown(False)
          self.do_SpinDown = 0

    def SpinUp_clicked(self):
        if ( self.do_SpinUp == 0):
          self.pushButtonSpinUp.setDown(True)
          self.do_SpinUp = 1
        else:
          self.pushButtonSpinUp.setDown(False)
          self.do_SpinUp = 0


    def S_clicked(self):
        if ( self.do_S == 0):
          self.DisableAllLandM()
          self.pushButtonS.setDown(True)
          self.do_S = 1
          self.pushButtonS0.setDown(True)
          self.do_S0 = 1
        else:

          self.pushButtonS.setDown(False)
          self.do_S = 0
          self.pushButtonS0.setDown(False)
          self.do_S0 = 0

    def P_clicked(self):
        if ( self.do_P == 0):
          self.DisableAllLandM()

          self.pushButtonP.setDown(True)
          self.do_P = 1

          self.pushButtonP_1.setDown(True)
          self.do_P_1 = 1
          self.pushButtonP0.setDown(True)
          self.do_P0 = 1
          self.pushButtonP1.setDown(True)
          self.do_P1 = 1
        else:

          self.pushButtonP.setDown(False)
          self.do_P = 0

          self.pushButtonP_1.setDown(False)
          self.do_P_1 = 0
          self.pushButtonP0.setDown(False)
          self.do_P0 = 0
          self.pushButtonP1.setDown(False)
          self.do_P1 = 0


    def D_clicked(self):
        if ( self.do_D == 0):
          self.DisableAllLandM()

          self.pushButtonD.setDown(True)
          self.do_D = 1

          self.pushButtonD_2.setDown(True)
          self.do_D_2 = 1
          self.pushButtonD_1.setDown(True)
          self.do_D_1 = 1
          self.pushButtonD0.setDown(True)
          self.do_D0 = 1
          self.pushButtonD1.setDown(True)
          self.do_D1 = 1
          self.pushButtonD2.setDown(True)
          self.do_D2 = 1

        else:

          self.pushButtonD.setDown(False)
          self.do_D = 0

          self.pushButtonD_2.setDown(False)
          self.do_D_2 = 0
          self.pushButtonD_1.setDown(False)
          self.do_D_1 = 0
          self.pushButtonD0.setDown(False)
          self.do_D0 = 0
          self.pushButtonD1.setDown(False)
          self.do_D1 = 0
          self.pushButtonD2.setDown(False)
          self.do_D2 = 0

    def F_clicked(self):
        if ( self.do_F == 0):
          self.DisableAllLandM()

          self.pushButtonF.setDown(True)
          self.do_F = 1

          self.pushButtonF_3.setDown(True)
          self.do_F_3 = 1
          self.pushButtonF_2.setDown(True)
          self.do_F_2 = 1
          self.pushButtonF_1.setDown(True)
          self.do_F_1 = 1
          self.pushButtonF0.setDown(True)
          self.do_F0 = 1
          self.pushButtonF1.setDown(True)
          self.do_F1 = 1
          self.pushButtonF2.setDown(True)
          self.do_F2 = 1
          self.pushButtonF3.setDown(True)
          self.do_F3 = 1

        else:

          self.pushButtonF.setDown(False)
          self.do_F = 0

          self.pushButtonF_3.setDown(False)
          self.do_F_3 = 0
          self.pushButtonF_2.setDown(False)
          self.do_F_2 = 0
          self.pushButtonF_1.setDown(False)
          self.do_F_1 = 0
          self.pushButtonF0.setDown(False)
          self.do_F0 = 0
          self.pushButtonF1.setDown(False)
          self.do_F1 = 0
          self.pushButtonF2.setDown(False)
          self.do_F2 = 0
          self.pushButtonF3.setDown(False)
          self.do_F3 = 0

# m buttons

    def S0_clicked(self):
        if ( self.do_S0 == 0):

          self.pushButtonS0.setDown(True)
          self.do_S0 = 1

          self.pushButtonS.setDown(True)

        else:

          self.pushButtonS0.setDown(False)
          self.do_S0 = 0


    def P_1_clicked(self):
        if ( self.do_P_1 == 0):

          self.pushButtonP_1.setDown(True)
          self.do_P_1 = 1

          self.pushButtonP.setDown(True)

        else:

          self.pushButtonP_1.setDown(False)
          self.do_P_1 = 0

          if(self.do_P_1==0 and self.do_P0==0 and self.do_P1==0):
            self.pushButtonP.setDown(False)

    def P0_clicked(self):
        if ( self.do_P0 == 0):

          self.pushButtonP0.setDown(True)
          self.do_P0 = 1

          self.pushButtonP.setDown(True)

        else:

          self.pushButtonP0.setDown(False)
          self.do_P0 = 0

          if(self.do_P_1==0 and self.do_P0==0 and self.do_P1==0):
            self.pushButtonP.setDown(False)

    def P1_clicked(self):
        if ( self.do_P1 == 0):

          self.pushButtonP1.setDown(True)
          self.do_P1 = 1

          self.pushButtonP.setDown(True)

        else:

          self.pushButtonP1.setDown(False)
          self.do_P1 = 0

          if(self.do_P_1==0 and self.do_P0==0 and self.do_P1==0):
            self.pushButtonP.setDown(False)


    def D_2_clicked(self):
        if ( self.do_D_2 == 0):

          self.pushButtonD_2.setDown(True)
          self.do_D_2 = 1

          self.pushButtonD.setDown(True)

        else:

          self.pushButtonD_2.setDown(False)
          self.do_D_2 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setDown(False)

    def D_1_clicked(self):
        if ( self.do_D_1 == 0):

          self.pushButtonD_1.setDown(True)
          self.do_D_1 = 1

          self.pushButtonD.setDown(True)

        else:

          self.pushButtonD_1.setDown(False)
          self.do_D_1 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setDown(False)

    def D0_clicked(self):
        if ( self.do_D0 == 0):

          self.pushButtonD0.setDown(True)
          self.do_D0 = 1

          self.pushButtonD.setDown(True)

        else:

          self.pushButtonD0.setDown(False)
          self.do_D0 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setDown(False)

    def D1_clicked(self):
        if ( self.do_D1 == 0):

          self.pushButtonD1.setDown(True)
          self.do_D1 = 1

          self.pushButtonD.setDown(True)

        else:

          self.pushButtonD1.setDown(False)
          self.do_D1 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setDown(False)

    def D2_clicked(self):
        if ( self.do_D2 == 0):

          self.pushButtonD2.setDown(True)
          self.do_D2 = 1

          self.pushButtonD.setDown(True)


        else:

          self.pushButtonD2.setDown(False)
          self.do_D2 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setDown(False)


    def F_3_clicked(self):
        if ( self.do_F_3 == 0):

          self.pushButtonF_3.setDown(True)
          self.do_F_3 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF_3.setDown(False)
          self.do_F_3 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)


    def F_2_clicked(self):
        if ( self.do_F_2 == 0):

          self.pushButtonF_2.setDown(True)
          self.do_F_2 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF_2.setDown(False)
          self.do_F_2 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)


    def F_1_clicked(self):
        if ( self.do_F_1 == 0):

          self.pushButtonF_1.setDown(True)
          self.do_F_1 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF_1.setDown(False)
          self.do_F_1 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)


    def F0_clicked(self):
        if ( self.do_F0 == 0):

          self.pushButtonF0.setDown(True)
          self.do_F0 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF0.setDown(False)
          self.do_F0 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)


    def F1_clicked(self):
        if ( self.do_F1 == 0):

          self.pushButtonF1.setDown(True)
          self.do_F1 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF1.setDown(False)
          self.do_F1 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)


    def F2_clicked(self):
        if ( self.do_F2 == 0):

          self.pushButtonF2.setDown(True)
          self.do_F2 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF2.setDown(False)
          self.do_F2 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)


    def F3_clicked(self):
        if ( self.do_F3 == 0):

          self.pushButtonF3.setDown(True)
          self.do_F3 = 1

          self.pushButtonF.setDown(True)

        else:

          self.pushButtonF3.setDown(False)
          self.do_F3 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setDown(False)



    def printButtonPressed(self):
        # This is executed when the button is pressed
        print('printButtonPressed')

app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()
