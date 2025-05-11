# QtBandViewer
# View Band Structures with unfolding and atomic orbital character
# using PyQt and MatplotLib.
# Version 0.6
# Copyright Carlos Loia Reis /INESC-MN 2021.

import sys, os

from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QHBoxLayout, QGridLayout, QPushButton
from PyQt5.QtWidgets import QMessageBox, QDialog,QTextEdit

#JLM 10/05/2025

import matplotlib
matplotlib.use('QtAgg')


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from BandDat import BandDat_t, DosDat_t


class OrbInfoDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        # Load the dialog's GUI
        uic.loadUi("orbinfo.ui", self)


class Ui(QtWidgets.QMainWindow):

    def select_iselected(self,iband):
      self.current_alpha[:] = 0.0
      for iorb in range(self.bdat.nbaslcao):
        if (self.iselect[iorb] == 1) :
#          print ("selecting", iorb)
          self.current_alpha[:]=self.current_alpha[:]+ self.bdat.basxpsi[:,iband,iorb]

    def select_l(self,l,iband):
      self.current_alpha[:] = 0.0
      for iorb in range(self.bdat.nbaslcao):
        if self.bdat.infolcao[iorb,3] == l :
          self.current_alpha[:]=self.current_alpha[:]+ self.bdat.basxpsi[:,iband,iorb]

    def select_species(self,ispecies,iband):
      self.current_alpha[:] = 0.0
      for iorb in range(self.bdat.nbaslcao):
        if self.bdat.infolcao[iorb,0] == ispecies :
          self.current_alpha[:]=self.current_alpha[:]+ self.bdat.basxpsi[:,iband,iorb]


    def plt_update(self):

#      print("plt_update irk min max", self.irk_min, self.irk_max)

      self.point_size =self.spinBoxPointSize.value()

      for i in range(self.bdat.neig):
        if (self.removed[i] == False) :
#          print("removing ", i)
          self.active_plots[i].remove()
          self.plots[i].remove()

          self.removed[i]=True

      for i in range(self.iband_min,self.iband_max):
        self.removed[i]=False
        if(self.select_mode ==0):
#          print ("selecting l=",self.selected_l)
          self.select_l(self.selected_l,i)

        if(self.select_mode ==1):
#          print ("selecting species=",self.selected_species)
          self.select_species(self.selected_species,i)

        if(self.select_mode ==-1):
          self.select_iselected(i)

        self.active_plots[i] = self.ax1.scatter(
        self.bdat.xk[self.irk_min:self.irk_max],
        self.bdat.bands[self.irk_min:self.irk_max,i],
        alpha=self.current_alpha[self.irk_min:self.irk_max]*self.bdat.pkn[i,self.irk_min:self.irk_max],
        s = self.point_size,
        color=self.orb_color[self.irk_min:self.irk_max,i,:],
        edgecolor='none')

        if (self.ev_max[i] > 0.0) :

          self.plots[i] = self.ax1.scatter(
          self.bdat.xk[self.irk_min:self.irk_max],
          self.bdat.bands[self.irk_min:self.irk_max,i],
          alpha=self.bdat.pkn[i,self.irk_min:self.irk_max],
          s = self.point_size,
          color="Red",
          edgecolor='none',
          visible = self.plot_visible )

        else:

          self.plots[i] = self.ax1.scatter(
          self.bdat.xk[self.irk_min:self.irk_max],
          self.bdat.bands[self.irk_min:self.irk_max,i],
          alpha=self.bdat.pkn[i,self.irk_min:self.irk_max],
          s = self.point_size,
          color="Blue",
          edgecolor='none',
          visible = self.plot_visible )

      if not self.has_dos :
        plt.draw()
        return

      self.current_dos[:] = 0.0
      for iorb in range(self.bdat.nbaslcao):
        if (self.iselect[iorb] == 1) :
          self.current_dos[:] = self.current_dos[:] + self.ddat.projected_dos[iorb,:]


      self.ax2.clear()
      self.ax2.set_xlabel('DOS [el/cell/eV]')
      self.ax2.yaxis.set_ticklabels([])
      self.ax2.grid(True)

      self.ax2.set_xlim(0.0,max(self.ddat.total_dos)*1.05)
      (ymin,ymax) = self.ax1.get_ylim()
      self.ax2.set_ylim(ymin,ymax)
      self.ax2.locator_params(axis="x", nbins=4)

      self.dos_plot=self.ax2.plot(self.current_dos,self.ddat.energy,color="Red")
      self.ax2.fill_betweenx(self.ddat.energy,0.0000,self.current_dos, facecolor='Red', alpha=0.5)


      self.all_dos_plot=self.ax2.plot(self.ddat.total_dos,self.ddat.energy,color="Grey")
      self.ax2.fill_betweenx(self.ddat.energy,0.0000,self.ddat.total_dos, facecolor='Grey', alpha=0.5)
      plt.draw()
      return


    def on_xlims_change(self,event_ax):
      print("updated xlims: ", event_ax.get_xlim())
      (xmin, xmax) = event_ax.get_xlim()
      print ("xmin", xmin)
      print ("xmax", xmax)
      self.irk_min = 0
      self.irk_max = self.bdat.nrk

      for irk in range(self.bdat.nrk):
        if (self.bdat.xk[irk] > xmin):
          self.irk_min = irk
          break

      for irk in range(self.bdat.nrk):
        if (self.bdat.xk[irk] > xmax):
          self.irk_max = irk
          self.irk_max = min(irk,self.bdat.nrk)
          break

#      self.irk_max = self.irk_max -1

      print ("irk min", self.irk_min)
      print ("xk min", self.bdat.xk[self.irk_min])
      print ("irk max",self.irk_max)
      print ("xk max", self.bdat.xk[self.irk_max-1])
      self.plt_update()


    def SetBandExtrema(self,ymin,ymax):
      self.iband_min = 0
      self.iband_max = self.bdat.neig

      for i in range(self.bdat.neig):
        self.ev_max[i] = max(self.bdat.bands[self.irk_min:self.irk_max,i])
        self.ev_min[i] = min(self.bdat.bands[self.irk_min:self.irk_max,i])


      for iband in range(self.bdat.neig):
        if (self.ev_min[iband] < ymin):
          self.iband_min = iband
#          break

      for iband in range(self.bdat.neig):
        if (self.ev_max[iband] < ymax):
          self.iband_max = iband
#          break

      self.iband_max= self.iband_max+1

#      print("iband_min, iband_max", self.iband_min,self.iband_max)

      self.iband_max = min(self.iband_max+6,self.bdat.neig)
      self.iband_min = max(self.iband_min-6,0)

#      self.iband_max = min(self.iband_max,self.bdat.neig)
#     self.iband_min = max(self.iband_min,0)


#      self.iband_max = self.bdat.neig
#      self.iband_min = 0



    def on_ylims_change(self,event_ax):
      print("updated ylims: ", event_ax.get_ylim())
      (ymin, ymax) = event_ax.get_ylim()
      self.SetBandExtrema(ymin,ymax)
      if self.has_dos:
        self.ax2.set_ylim(ymin,ymax)

      self.plt_update()


    def on_ylims_change_old(self,event_ax):
#      print("updated ylims: ", event_ax.get_ylim())
      (ymin, ymax) = event_ax.get_ylim()
      self.iband_min = 0
      self.iband_max = self.bdat.neig

      for i in range(self.bdat.neig):
        self.ev_max[i] = max(self.bdat.bands[self.irk_min:self.irk_max,i])
        self.ev_min[i] = min(self.bdat.bands[self.irk_min:self.irk_max,i])


      for iband in range(self.bdat.neig):
        if (self.ev_min[iband] > ymin):
          self.iband_min = iband
          break

      for iband in range(self.bdat.neig):
        if (self.ev_max[iband]> ymax):
          self.iband_max = iband
          break

      self.iband_max= self.iband_max+1

#      print("iband_min, iband_max", self.iband_min,self.iband_max)
      self.plt_update()


    def toggleBands(self):
      if self.toggle_mode == 0 :
         self.toggle_mode =1
         self.pushButtonToggle.setText("All Bands")
      else:
         self.toggle_mode =0
         self.pushButtonToggle.setText("Proj. Bands")


#      for i in range(self.bdat.neig):

      for i in range(self.iband_min,self.iband_max):
        self.active_plots[i].set_visible(not self.active_plots[i].get_visible() )
        self.plots[i].set_visible(not self.plots[i].get_visible() )
        self.plot_visible = self.plots[i].get_visible()

        if(self.select_mode ==0):
#          print ("selecting l=",self.selected_l)
          self.select_l(self.selected_l,i)

        if(self.select_mode ==1):
#          print ("selecting species=",self.selected_species)
          self.select_species(self.selected_species,i)

        if(self.select_mode ==-1):
          self.select_iselected(i)

      if self.has_dos :
        self.current_dos[:] = 0.0
        for iorb in range(self.bdat.nbaslcao):
          if (self.iselect[iorb] == 1) :
            self.current_dos[:] = self.current_dos[:] + self.ddat.projected_dos[iorb,:]

      plt.draw()
      return


    def toggle_images(self,event):
#      print("Hello ", self.bdat.title)
      if event.key =='l':
#        print ("select_mode=angular momentum")
        self.select_mode=0
        return
      if event.key =='z':
#        print ("select_mode=atomic species")
        self.select_mode=1
        return
      if event.key =='s':
#        print ("selected_l=0")
        self.selected_l =0
        self.plt_update()
        return
      if event.key =='p':
#        print ("selected_l=1")
        self.selected_l = 1
        self.plt_update()
        return
      if event.key =='d':
#        print ("selected_l=2")
        self.selected_l = 2
        self.plt_update()
        return
      if event.key =='c':
        self.orb_color_mode = not self.orb_color_mode
        self.plt_update()
        return
      if event.key =='1':
#        print ("selected_species=1")
        self.selected_species =1
        self.plt_update()
        return
      if event.key =='2':
#        print ("selected_species=2")
        self.selected_species =2
        self.plt_update()
        return
      if event.key =='3':
#        print ("selected_species=3")
        self.selected_species =3
        self.plt_update()
        return
      if event.key =='4':
#        print ("selected_species=4")
        self.selected_species =4
        self.plt_update()
        return

      if event.key == 't':
#        print("hahah")
#        print("select_mode selected", self.select_mode, self.selected_l)
         self.toggleBands()
         return

    def close(self):
      print ("closing")
      super(Ui, self).close()

    def __init__(self):

        self.cwd = os.getcwd();

        super(Ui, self).__init__()
        uic.loadUi('QtBandViewer.ui', self)

        self.setGeometry(900,0,self.width(),self.height())

        self.show()

        self.groupBoxL.setEnabled(False)
        self.groupBoxM.setEnabled(False)
        self.groupBoxSO.setEnabled(False)
        self.spinBoxPointSize.setEnabled(False)
        self.groupBoxAtoms.setEnabled(False)
        self.groupBoxSpecies.setEnabled(False)
        self.pushButtonFullRange.setEnabled(False)
        self.pushButtonS.setEnabled(False)
        self.pushButtonS0.setEnabled(False)


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

#global select_mode
        self.select_mode=-1
        self.selected_l = 0
        self.selected_species = 1
        self.toggle_mode=0
        self.orb_color_mode=True

        self.checkBoxAtoms =[None]*80
        self.checkBoxAtoms[0] = self.checkBoxAtom1
        self.checkBoxAtoms[0] = self.checkBoxAtom1
        self.checkBoxAtoms[1] = self.checkBoxAtom2
        self.checkBoxAtoms[2] = self.checkBoxAtom3
        self.checkBoxAtoms[3] = self.checkBoxAtom4
        self.checkBoxAtoms[4] = self.checkBoxAtom5
        self.checkBoxAtoms[5] = self.checkBoxAtom6
        self.checkBoxAtoms[6] = self.checkBoxAtom7
        self.checkBoxAtoms[7] = self.checkBoxAtom8
        self.checkBoxAtoms[8] = self.checkBoxAtom9
        self.checkBoxAtoms[9] = self.checkBoxAtom10
        self.checkBoxAtoms[10] = self.checkBoxAtom11
        self.checkBoxAtoms[11] = self.checkBoxAtom12
        self.checkBoxAtoms[12] = self.checkBoxAtom13
        self.checkBoxAtoms[13] = self.checkBoxAtom14
        self.checkBoxAtoms[14] = self.checkBoxAtom15
        self.checkBoxAtoms[15] = self.checkBoxAtom16
        self.checkBoxAtoms[16] = self.checkBoxAtom17
        self.checkBoxAtoms[17] = self.checkBoxAtom18
        self.checkBoxAtoms[18] = self.checkBoxAtom19
        self.checkBoxAtoms[19] = self.checkBoxAtom20
        self.checkBoxAtoms[20] = self.checkBoxAtom21
        self.checkBoxAtoms[21] = self.checkBoxAtom22
        self.checkBoxAtoms[22] = self.checkBoxAtom23
        self.checkBoxAtoms[23] = self.checkBoxAtom24
        self.checkBoxAtoms[24] = self.checkBoxAtom25
        self.checkBoxAtoms[25] = self.checkBoxAtom26
        self.checkBoxAtoms[26] = self.checkBoxAtom27
        self.checkBoxAtoms[27] = self.checkBoxAtom28
        self.checkBoxAtoms[28] = self.checkBoxAtom29
        self.checkBoxAtoms[29] = self.checkBoxAtom30
        self.checkBoxAtoms[30] = self.checkBoxAtom31
        self.checkBoxAtoms[31] = self.checkBoxAtom32
        self.checkBoxAtoms[32] = self.checkBoxAtom33
        self.checkBoxAtoms[33] = self.checkBoxAtom34
        self.checkBoxAtoms[34] = self.checkBoxAtom35
        self.checkBoxAtoms[35] = self.checkBoxAtom36
        self.checkBoxAtoms[36] = self.checkBoxAtom37
        self.checkBoxAtoms[37] = self.checkBoxAtom38
        self.checkBoxAtoms[38] = self.checkBoxAtom39
        self.checkBoxAtoms[39] = self.checkBoxAtom40
        self.checkBoxAtoms[40] = self.checkBoxAtom41
        self.checkBoxAtoms[41] = self.checkBoxAtom42
        self.checkBoxAtoms[42] = self.checkBoxAtom43
        self.checkBoxAtoms[43] = self.checkBoxAtom44
        self.checkBoxAtoms[44] = self.checkBoxAtom45
        self.checkBoxAtoms[45] = self.checkBoxAtom46
        self.checkBoxAtoms[46] = self.checkBoxAtom47
        self.checkBoxAtoms[47] = self.checkBoxAtom48
        self.checkBoxAtoms[48] = self.checkBoxAtom49
        self.checkBoxAtoms[49] = self.checkBoxAtom50
        self.checkBoxAtoms[50] = self.checkBoxAtom51
        self.checkBoxAtoms[51] = self.checkBoxAtom52
        self.checkBoxAtoms[52] = self.checkBoxAtom53
        self.checkBoxAtoms[53] = self.checkBoxAtom54
        self.checkBoxAtoms[54] = self.checkBoxAtom55
        self.checkBoxAtoms[55] = self.checkBoxAtom56
        self.checkBoxAtoms[56] = self.checkBoxAtom57
        self.checkBoxAtoms[57] = self.checkBoxAtom58
        self.checkBoxAtoms[58] = self.checkBoxAtom59
        self.checkBoxAtoms[59] = self.checkBoxAtom60
        self.checkBoxAtoms[60] = self.checkBoxAtom61
        self.checkBoxAtoms[61] = self.checkBoxAtom62
        self.checkBoxAtoms[62] = self.checkBoxAtom63
        self.checkBoxAtoms[63] = self.checkBoxAtom64
        self.checkBoxAtoms[64] = self.checkBoxAtom65
        self.checkBoxAtoms[65] = self.checkBoxAtom66
        self.checkBoxAtoms[66] = self.checkBoxAtom67
        self.checkBoxAtoms[67] = self.checkBoxAtom68
        self.checkBoxAtoms[68] = self.checkBoxAtom69
        self.checkBoxAtoms[69] = self.checkBoxAtom70
        self.checkBoxAtoms[70] = self.checkBoxAtom71
        self.checkBoxAtoms[71] = self.checkBoxAtom72
        self.checkBoxAtoms[72] = self.checkBoxAtom73
        self.checkBoxAtoms[73] = self.checkBoxAtom74
        self.checkBoxAtoms[74] = self.checkBoxAtom75
        self.checkBoxAtoms[75] = self.checkBoxAtom76
        self.checkBoxAtoms[76] = self.checkBoxAtom77
        self.checkBoxAtoms[77] = self.checkBoxAtom78
        self.checkBoxAtoms[78] = self.checkBoxAtom79
        self.checkBoxAtoms[79] = self.checkBoxAtom80

        for i in range(80) :
          self.checkBoxAtoms[i].hide()

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

      self.Run()

    def InitList(self):

        if (self.lso) :
          self.pushButtonSpinUp.setDown(True)
          self.pushButtonSpinDown.setDown(True)
        else:
          self.pushButtonSpinUp.setDown(False)
          self.pushButtonSpinDown.setDown(False)

        self.info_itype = [None]*self.bdat.nbaslcao
        self.info_iatom = [None]*self.bdat.nbaslcao
        self.info_n     = [None]*self.bdat.nbaslcao
        self.info_l     = [None]*self.bdat.nbaslcao
        self.info_m     = [None]*self.bdat.nbaslcao

        for iorb in range(0,self.bdat.nbaslcao):
          self.info_itype[iorb] = self.bdat.infolcao[iorb,0]
          self.info_iatom[iorb]=self.bdat.infolcao[iorb,1]
          self.info_n[iorb]=self.bdat.infolcao[iorb,2]
          self.info_l[iorb]=self.bdat.infolcao[iorb,3]
          self.info_m[iorb]=self.bdat.infolcao[iorb,4]

        self.nspecies = self.bdat.ntype

        self.nameat =[None]*(self.nspecies)

        self.name_keys = {}

        for ispecies in range(0,self.nspecies):
          self.nameat[ispecies] = self.bdat.nameat[ispecies]
          self.name_keys[self.nameat[ispecies]] = ispecies



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



#        print("self.nspecies is ", self.nspecies)

#        print("self.nameat is ", self.nameat)

        self.natom=[None]*(self.nspecies)


        for ispecies in range(0,self.nspecies):
          self.natom[ispecies]=0
#          print("ispecies self.natom", ispecies,self.natom)
          if (ispecies ==0):
            self.lineEditSpecies.setText( str(ispecies+1))
          if (ispecies >0):
            self.lineEditSpecies.setText( self.lineEditSpecies.text() +"," + str(ispecies+1))

        for iorb in range(0,self.bdat.nbaslcao):
          ispecies = self.info_itype[iorb] - 1
          if (self.info_iatom[iorb] > self.natom[ispecies]) :
            self.natom[ispecies]=self.natom[ispecies] +1


        self.atom_list_str = "";
        iatom_cb=0

        for ispecies in range(0, self.nspecies):
          for iatom in range(0, self.natom[ispecies]):
#            print ("iatom=", iatom)
            self.atom_list_str = self.atom_list_str + str(iatom+1) + " "
            if (iatom +1 < self.natom[ispecies]):
              self.atom_list_str = self.atom_list_str + ","
            if (iatom +1 == self.natom[ispecies]):
              self.atom_list_str =self.atom_list_str + "\n"


#            print(self.nameat[ispecies])
            cb_txt = self.nameat[ispecies] + str(iatom+1)
            self.checkBoxAtoms[iatom_cb].show()
            self.checkBoxAtoms[iatom_cb].setChecked(True)
            self.checkBoxAtoms[iatom_cb].setText(cb_txt)

            iatom_cb = iatom_cb+1


#        self.textEditAtoms.setPlainText(self.atom_list_str)

#        print("self.natom list", self.natom)

        self.lineEditSpecies.setEnabled(False)
#        self.textEditAtoms.setEnabled(False)


    def CheckSpecies(self,icheck):
      if (icheck ==0 ) :
        self.lineEditSpecies.setEnabled(True)
      else:
        self.lineEditSpecies.setEnabled(False)


    def Load(self):

        fname_info = QFileDialog.getOpenFileName(self,'Open file',".","Band info binary files  (BAND*DAT*bv)" )

        dirname =  os.path.dirname(fname_info[0]) ## directory of file

        print("dirname",dirname)

        # dos stuff
        self.has_dos = False
        if os.path.isfile(dirname+"/dos.bv"):
          print ("found dos.bv")
          self.has_dos = True
        else:
          print ("no information about DOS")

        if self.has_dos:
          gs = gridspec.GridSpec(1, 2,width_ratios=[4,1])
          self.ax1 = plt.subplot(gs[0])
          self.ax2 = plt.subplot(gs[1])
        else:
          gs = gridspec.GridSpec(1, 1)
          self.ax1 = plt.subplot(gs[0])

        self.annot = self.ax1.annotate("", xy=(0,0), xytext=(5,5),textcoords="figure pixels",
                    #bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))

        self.annot.set_visible(False)

        self.point_size =self.spinBoxPointSize.value()

        if self.has_dos:

          self.ddat = DosDat_t()

          self.ddat.ReadFile(dirname+"/dos.bv")
          self.current_dos  = np.empty([self.ddat.n],dtype=np.float64)


          self.dos_active_plots=[None]*self.ddat.norb
          self.dos_plot        =[None]
          self.dos_removed     =[bool]*self.ddat.norb


          self.ax2.set_xlim(0.0,max(self.ddat.total_dos)*1.05)
          self.ax2.locator_params(axis="x", nbins=4)
#          self.ax2.get_yaxis().set_visible(False)
          self.ax2.yaxis.set_ticklabels([])
#          self.ax2.grid(True)
          self.ax2.set_xlabel('DOS [el/cell/eV]')
          self.all_dos_plot=self.ax2.plot(self.ddat.total_dos,self.ddat.energy,color="Grey")
          self.ax2.fill_betweenx(self.ddat.energy,0.0000,self.ddat.total_dos, facecolor='Grey', alpha=0.5)



        band_lines_name = dirname + "/BAND_LINES.DAT"

        print("band_lines_name",band_lines_name)

        fin_band_lines = open(band_lines_name, "r")

        line = fin_band_lines.readline().split()

        self.npanel = int(line[0])
        nbands      = int(line[1])

        print(self.npanel)
        print(nbands)

        self.kvec_start = np.zeros((3, self.npanel), dtype=np.float32)
        self.kvec_end = np.zeros((3, self.npanel), dtype=np.float32)

        self.klabel_start =[]*self.npanel
        self.klabel_end   =[]*self.npanel

        for i in range(self.npanel) :
          spl =fin_band_lines.readline().split()
          self.kvec_start[0,i] = float(spl[0])
          self.kvec_start[1,i] = float(spl[1])
          self.kvec_start[2,i] = float(spl[2])
          self.kvec_end[0,i]   = float(spl[3])
          self.kvec_end[1,i]   = float(spl[4])
          self.kvec_end[2,i]   = float(spl[5])
          self.klabel_start.append(spl[7].replace("Gamma","G"))
          self.klabel_end.append(spl[9].replace("Gamma","G"))

          print("strt:", self.kvec_start[:,i])
          print("end:",  self.kvec_end[:,i])

          print("strt:", self.klabel_start[i])
          print("end:",  self.klabel_end[i])


        self.kbutton = []*self.npanel

        self.hbox_kpanel = QGridLayout()

        self.groupBoxKPanel.setLayout(self.hbox_kpanel)

        self.kbutton_all  = QPushButton("All")
        self.kbutton_all.clicked.connect(self.kbutton_all_clicked)

        self.kbutton_none = QPushButton("None")
        self.kbutton_none.clicked.connect(self.kbutton_none_clicked)

        iy=0
        ix = 0
        self.hbox_kpanel.addWidget(self.kbutton_all,iy,ix)

        iy=0
        ix = 1

        self.hbox_kpanel.addWidget(self.kbutton_none,iy,ix)

        ix = 0
        iy = 1

        for i in range(self.npanel) :
          self.kbutton.append(QPushButton(self.klabel_start[i]+ "-" + self.klabel_end[i]))
          self.kbutton[i].clicked.connect(self.kpanelClicked)
          self.kbutton[i].setCheckable(True)
          self.kbutton[i].setChecked(True)

          self.hbox_kpanel.addWidget(self.kbutton[i],iy,ix)
          ix = ix+1

          print("i,ix,iy", i,ix,iy)

          if (ix == 5) :
            iy = iy+1
            ix = 0

        if fname_info[0] == "" : return

        self.fname = fname_info[0]

        so_str = "SO"

        if so_str in self.fname:
          self.lso = True
#          print ("Found!")

        else:
          self.lso = False
#          print ("Not found!")

        if self.lso :
          self.groupBoxSO.setEnabled(True)
        else:
          self.groupBoxSO.setEnabled(False)


#        quit()

      ## Begin main code.

        self.bdat = BandDat_t();
        self.bdat.ReadFileFl32(self.fname)

        # print infolcao
#!<  information about the original atomic orbital.
#(type of atom, atom of that type, n,l,m)

        for iorb in range(self.bdat.nbaslcao):
          print("iorb %3d itype %3d iatom %d l=%2d m=%2d "
          %(iorb, self.bdat.infolcao[iorb,0],
          self.bdat.infolcao[iorb,1],
          self.bdat.infolcao[iorb,3],
          self.bdat.infolcao[iorb,4]))

        self.irk_min=0
        self.irk_max=self.bdat.nrk

        self.delta_hx = (self.bdat.xmax-self.bdat.xmin)/self.bdat.nrk/2


        self.ev_min = np.zeros([self.bdat.neig])
        self.ev_max = np.zeros([self.bdat.neig])

        self.iband_min = 0
        self.iband_max = self.bdat.neig

        for i in range(self.iband_min,self.iband_max):
          self.ev_max[i] = max(self.bdat.bands[:,i])
          self.ev_min[i] = min(self.bdat.bands[:,i])

        self.default_ymin = -1.0
        self.default_ymax = +2.0

        self.ax1.set_xlim(self.bdat.xmin,self.bdat.xmax)
        self.ax1.set_ylim(self.default_ymin,self.default_ymax)
        if self.has_dos:
          self.ax2.set_ylim(self.default_ymin,self.default_ymax)
          self.ax2.grid(True)

        self.SetBandExtrema(self.default_ymin,self.default_ymax)

#        plt.xticks(self.bdat.xklab,self.bdat.labels)

        self.ax1.set_xticks(self.bdat.xklab)
        self.ax1.set_xticklabels(self.bdat.labels)
#        plt.grid(True)
#        plt.ylabel('E [eV]')
        self.ax1.grid(True)
        self.ax1.set_ylabel('E [eV]')

#        plt.title(self.bdat.subtitle.decode("UTF-8").lstrip().rstrip())
        plt.suptitle(self.bdat.title.decode("UTF-8").lstrip().rstrip())
        self.ax1.set_title(self.bdat.subtitle.decode("UTF-8").lstrip().rstrip())
#        plt.legend()

# Test code
        self.orb_color = np.zeros([self.bdat.nrk,self.bdat.neig,3])
        for iorb in range(self.bdat.nbaslcao):
          if (self.bdat.infolcao[iorb][3] == 0) :
            self.orb_color[:,:,0] = self.orb_color[:,:,0] +  self.bdat.basxpsi[:,:,iorb]
          if (self.bdat.infolcao[iorb][3] == 1) :
            self.orb_color[:,:,2] = self.orb_color[:,:,2] +  self.bdat.basxpsi[:,:,iorb]
          if (self.bdat.infolcao[iorb][3] == 2) :
            self.orb_color[:,:,1] = self.orb_color[:,:,1] +  self.bdat.basxpsi[:,:,iorb]


        self.point_size =self.spinBoxPointSize.value()

        self.active_plots=[None]*self.bdat.neig
        self.plots       =[None]*self.bdat.neig
        self.removed =[bool]*self.bdat.neig

        for i in range(self.iband_min,self.iband_max):

          self.active_plots[i] = self.ax1.scatter(self.bdat.xk, self.bdat.bands[:,i],
                      alpha=self.bdat.pkn[i,:],
                      s = self.point_size,
                      color=self.orb_color[:,i,:],edgecolor='none')

          if (self.ev_max[i] > 0.0) :
            self.plots[i]        = self.ax1.scatter(self.bdat.xk, self.bdat.bands[:,i],
                      alpha=self.bdat.pkn[i,:],
                      s = self.point_size, color="Red",
                      edgecolor='none')
          else:
            self.plots[i]        = self.ax1.scatter(self.bdat.xk, self.bdat.bands[:,i],
                      alpha=self.bdat.pkn[i,:],
                      s = self.point_size, color="Blue",
                      edgecolor='none')


          self.active_plots[i].set_visible(False)


          self.removed[i] = False

        self.plot_visible = True
        self.current_alpha=np.zeros(self.bdat.nrk)

        plt.connect('key_press_event', self.toggle_images)
        plt.connect('motion_notify_event', self.hover)
        plt.connect('button_press_event', self.click)


#print(plt.rcParams)
#added JLM 10/05/2025.  The next four lines that control plot
#appearence are not commented in some versions
#but thay cause a crash...
#        plt.rcParams['keymap.save'].remove('s')
#        plt.rcParams['keymap.pan'].remove('p')
#        plt.rcParams['keymap.all_axes'].remove('a')
#        plt.rcParams['keymap.yscale'].remove('l')

        fig =plt.gcf()
        fig.set_size_inches(8,7)
        self.ax = plt.gca()
        plt.tight_layout()

#        self.plt_update()

        self.ax1.callbacks.connect('xlim_changed', self.on_xlims_change)
        self.ax1.callbacks.connect('ylim_changed', self.on_ylims_change)


        plt.show()


       ## end main code


        self.InitList()

        self.pushButtonRun.setEnabled(True)
        self.pushButtonToggle.setEnabled(True)

        self.iselect=[None]*self.bdat.nbaslcao

        self.groupBoxL.setEnabled(True)
        self.groupBoxM.setEnabled(True)
#        self.groupBoxSO.setEnabled(True)
        self.spinBoxPointSize.setEnabled(True)
        self.groupBoxAtoms.setEnabled(True)
        self.groupBoxSpecies.setEnabled(True)
        self.pushButtonFullRange.setEnabled(True)
        self.pushButtonS.setEnabled(True)
        self.pushButtonS0.setEnabled(True)


    def click(self,event):
      if not event.dblclick:
        return
      print (self.orb_info)
      info_dialog  = OrbInfoDialog(self)
      info_dialog.textEdit.setText(self.orb_info)
      info_dialog.exec()

    def hover(self,event):
      delta = self.delta_hx
      self.annot.set_visible(False)
#      plt.draw()
      if (event.xdata == None or event.ydata==None ) :
        return

#      first find nearest irk
      for irk in range(self.irk_min,self.irk_max):
        if (self.bdat.xk[irk] > event.xdata-delta
        and self.bdat.xk[irk] <event.xdata+delta ):
#          print ("irk is", irk,self.bdat.xk[irk], event.xdata)
          for ipanel in range(self.npanel):
            if (self.bdat.xk[irk] >= self.bdat.xk_start[ipanel]
            and self.bdat.xk[irk] <=  self.bdat.xk_end[ipanel]):
              print ("irk", irk, " is in panel #", ipanel)


          break

      delta = 0.01

      for iband in range(self.iband_min, self.iband_max):

        if (self.bdat.bands[irk,iband] > event.ydata -delta
        and self.bdat.bands[irk,iband] < event.ydata +delta):
#          print ("iband is", iband,self.bdat.bands[irk,iband], event.ydata)
          print("    itype | iatom |   l   |   m   |  proj  |")
          max_proj = -1


          self.orb_info = ""
          for iorb in range(self.bdat.nbaslcao):
            info_str =  "%7s %7d %7d %7d %10.3f\n"%(
            self.nameat[ self.bdat.infolcao[iorb,0] -1],
            self.bdat.infolcao[iorb,1],
            self.bdat.infolcao[iorb,3],
            self.bdat.infolcao[iorb,4],
            self.bdat.basxpsi[irk,iband,iorb])
            self.orb_info = self.orb_info + info_str

            if(self.bdat.basxpsi[irk,iband,iorb] > max_proj) :
              max_proj = self.bdat.basxpsi[irk,iband,iorb]
              iorb_max = iorb



          ev_info = "irk = %d\nkpt_fld   = (%8.5f,%8.5f,%8.5f)\nkpt_unfld = (%8.5f,%8.5f,%8.5f)\niband = %d, E = %8.5f (eV)"  % (
          irk,
          self.bdat.rkpt[irk,0],
          self.bdat.rkpt[irk,1],
          self.bdat.rkpt[irk,2],
          self.bdat.rkpt_fld[irk,0],
          self.bdat.rkpt_fld[irk,1],
          self.bdat.rkpt_fld[irk,2],
          iband,
          self.bdat.bands[irk,iband]
          )

          self.orb_info  = ev_info + "\n" + "    itype | iatom |   l   |   m   |  proj  |\n" + self.orb_info

          print (self.orb_info)


          print("maximum is iorb ",iorb_max)
          print("    itype | iatom |   l   |   m   |  proj  |")

#          prj_str = ""

          prj_str = "irk = %d, kpt_fld = (%8.5f,%8.5f,%8.5f), kpt_unfld = (%8.5f,%8.5f,%8.5f),\niband = %d, E = %8.5f (eV) %s # %d, l = %d, m = %d, proj=%6.3f" % (
          irk,
          self.bdat.rkpt[irk,0],
          self.bdat.rkpt[irk,1],
          self.bdat.rkpt[irk,2],
          self.bdat.rkpt_fld[irk,0],
          self.bdat.rkpt_fld[irk,1],
          self.bdat.rkpt_fld[irk,2],
          iband,
          self.bdat.bands[irk,iband],
          self.nameat[ self.bdat.infolcao[iorb_max,0] -1],
          self.bdat.infolcao[iorb_max,1],
          self.bdat.infolcao[iorb_max,3],
          self.bdat.infolcao[iorb_max,4],
          self.bdat.basxpsi[irk,iband,iorb_max])

          print (prj_str)
          self.annot.set_visible(True)
          self.annot.xycoords = "data"
          self.annot.textcoords = "offset points"
          print("JJJJJJJJJ")
          self.annot.xy = (event.xdata,event.ydata)
#          self.annot.xytext = (-20,200000)
          self.annot.set_text(prj_str)
          plt.draw()

#           break

    def FullRange(self):
#        plt.axis([self.bdat.xmin,self.bdat.xmax,self.bdat.ymin,self.bdat.ymax])
        self.ax1.set_xlim(self.bdat.xmin,self.bdat.xmax)
        self.ax1.set_ylim(self.bdat.ymin,self.bdat.ymax)
        if self.has_dos:
          self.ax2.set_ylim(self.bdat.ymin,self.bdat.ymax)

    def Run(self):

        if(self.toggle_mode ==0): self.toggleBands()

#      populate lineEditSpecies from check boxes

        self.atom_list_str = "";

#        self.textEditAtoms.setPlainText(self.atom_list_str)

        iatom_cb = 0
        for ispecies in range(0, self.nspecies):
          for iatom in range(0, self.natom[ispecies]):

            #print ("checked box ", ispecies,iatom,self.checkBoxAtoms[iatom_cb].isChecked())

            if self.checkBoxAtoms[iatom_cb].isChecked():
              self.atom_list_str = self.atom_list_str + str(iatom+1) + " "
              if (iatom +1 < self.natom[ispecies]):
                self.atom_list_str = self.atom_list_str + ","
              if (iatom +1 == self.natom[ispecies]):
                self.atom_list_str = self.atom_list_str + "\n"

            iatom_cb = iatom_cb+1

#        self.textEditAtoms.setPlainText(self.atom_list_str)


        for iorb in range(0,self.bdat.nbaslcao):
          self.iselect[iorb]=0


        SpeciesList  = self.lineEditSpecies.text().split(",");

#        print ("SpecList", SpeciesList, len(SpeciesList))

        if(SpeciesList[0]==""):
            return

#        print ("Species List is", SpeciesList)



#        AllAtomList  = self.textEditAtoms.toPlainText().split("\n");

        AllAtomList  = self.atom_list_str.split("\n");

#        print ("Atom List is", AllAtomList)


        for species in SpeciesList:
          ispecies = int(species)
          AtomList = AllAtomList[ispecies-1].split(",")
#          print ("iatom_list", AtomList )
          if not AtomList : print ("empty list?")
          for atom in AtomList:
#              print("atom in list", atom)
#            print ("atom",atom,"**")
            if atom=="" : break
            iatom = int(atom)
#              print("iatom", iatom)
            for iorb in range(0,self.bdat.nbaslcao):
              if(self.info_itype[iorb] == ispecies):
                if(self.info_iatom[iorb] == iatom):
                  if(self.info_l[iorb] == 0) :
                    if(self.do_S0 == 1):
                      self.iselect[iorb] = 1
                  if(self.info_l[iorb] == 1) :
                    if(self.info_m[iorb] == -1 and self.do_P_1 ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] ==  0 and self.do_P0  ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] ==  1 and self.do_P1  ==1):
                      self.iselect[iorb] = 1
                  if(self.info_l[iorb] == 2) :
                    if(self.info_m[iorb] == -2 and self.do_D_2 ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] == -1 and self.do_D_1 ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] == 0  and self.do_D0  ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] == 1  and self.do_D1  ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] == 2  and self.do_D2  ==1):
                      self.iselect[iorb] = 1
                  if(self.info_l[iorb] == 3) :
                    if(self.info_m[iorb] == -3 and self.do_F_3 ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] == -2 and self.do_F_2 ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] == -1 and self.do_F_1 ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] ==  0 and self.do_F0  ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] ==  1 and self.do_F1  ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] ==  2 and self.do_F2  ==1):
                      self.iselect[iorb] = 1
                    if(self.info_m[iorb] ==  3 and self.do_F3  ==1):
                      self.iselect[iorb] = 1


        if (self.lso) :
          for iorb in range(0,int(self.bdat.nbaslcao/2)):
#            print ("irob", iorb, self.do_SpinUp, self.do_SpinDown)
            if (self.do_SpinUp==0) :
              self.iselect[2*iorb]=0
            if (self.do_SpinDown==0) :
              self.iselect[2*iorb+1]=0

###debug
#        for iorb in range(0,self.bdat.nbaslcao):
#          print (iorb, "itype", self.info_itype[iorb], "iatom", self.info_iatom[iorb], "n", self.info_n[iorb], "l", self.info_l[iorb], "m", self.info_m[iorb], "self.iselect", self.iselect[iorb] )


        self.plt_update()

        return

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


#        print("ans ")
#        print(self.all_S)
#        print("ans ")
#        print(self.all_P)
#        print("ans ")
#        print(self.all_D)
#        print("ans ")
#        print(self.all_F)


    def DisableAllS(self):

        self.do_S = 0
        self.do_S0  = 0

        self.pushButtonS.setDown(False)
        self.pushButtonS0.setDown(False)

    def DisableAllP(self):

        self.do_P = 0

        self.do_P_1 = 0
        self.do_P0  = 0
        self.do_P1  = 0

        self.pushButtonP.setDown(False)
        self.pushButtonP_1.setDown(False)
        self.pushButtonP0.setDown(False)
        self.pushButtonP1.setDown(False)

    def DisableAllD(self):

        self.do_D = 0

        self.do_D_2  = 0
        self.do_D_1  = 0
        self.do_D0   = 0
        self.do_D1   = 0
        self.do_D2   = 0

        self.pushButtonD.setDown(False)
        self.pushButtonD_2.setDown(False)
        self.pushButtonD_1.setDown(False)
        self.pushButtonD0.setDown(False)
        self.pushButtonD1.setDown(False)
        self.pushButtonD2.setDown(False)


    def DisableAllF(self):

        self.do_F = 0

        self.do_F_3  = 0
        self.do_F_2  = 0
        self.do_F_1  = 0
        self.do_F0   = 0
        self.do_F1   = 0
        self.do_F2   = 0
        self.do_F3   = 0

        self.pushButtonF.setDown(False)
        self.pushButtonF_3.setDown(False)
        self.pushButtonF_2.setDown(False)
        self.pushButtonF_1.setDown(False)
        self.pushButtonF0.setDown(False)
        self.pushButtonF1.setDown(False)
        self.pushButtonF2.setDown(False)
        self.pushButtonF3.setDown(False)

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

        self.Run()

    def SpinUp_clicked(self):
        if ( self.do_SpinUp == 0):
          self.pushButtonSpinUp.setDown(True)
          self.do_SpinUp = 1
        else:
          self.pushButtonSpinUp.setDown(False)
          self.do_SpinUp = 0

        self.Run()


    def kbutton_all_clicked(self):
      for i in range(1,self.npanel) :
        self.kbutton[i].setChecked(True)

      xmin = self.bdat.xk_start[0]-0.0001
      xmax = self.bdat.xk_end[self.npanel-1]+0.0001
      self.ax1.set_xlim(xmin,xmax)

    def kbutton_none_clicked(self):
      print("none clicked")
      for i in range(1,self.npanel) :
        self.kbutton[i].setChecked(False)

      i=0
      self.kbutton[i].setChecked(True)
      xmin = self.bdat.xk_start[i]-0.0001
      xmax = self.bdat.xk_end[i]+0.0001
      self.ax1.set_xlim(xmin,xmax)

    def kpanelClicked(self):
      print("K panel was clicked");


      xmin = 1000
      xmax = -1

      icheck = 0

      for i in range (self.npanel):
        if (self.sender() == self.kbutton[i]) :
          print ("butons up to ", i," are pressed")
        if (self.kbutton[i].isChecked()):
          xmin = min(xmin,self.bdat.xk_start[i]-0.0001)
          xmax = max(xmax,self.bdat.xk_end[i]+0.0001)
          icheck = 1

      print("xmin xmax", xmin,xmax)
      if (icheck == 1) :
#        plt.xlim([xmin,xmax])
        self.ax1.set_xlim(xmin,xmax)
      else:
        xmin = 0.0
        xmax = 0.0
#        plt.xlim([xmin,xmax])

      for i in range (self.npanel):
        if (self.bdat.xk_end[i] < xmax and self.bdat.xk_start[i] > xmin) :
          self.kbutton[i].setChecked(True)


      # for i in range (j):
        # self.kbutton[i].setChecked(True)


      # for i in range (self.npanel):
        # if i <=k :
          # self.kbutton[i].setChecked(True)
        # else :
          # self.kbutton[i].setChecked(False)



      # for i in range (self.npanel):
        # if (self.sender() == self.kbutton[i]) :
          # self.kbutton[i].setChecked(True)
          # print("panel # ",i, " was clicked")
# #          plt.xlim([self.bdat.xk_start[i]-0.001,self.bdat.xk_end[i]+0.001])

        # else:
          # self.kbutton[i].setChecked(False)



    def S_clicked(self):
        if ( self.do_S == 0):
          self.DisableAllS()
#          self.pushButtonS.setChecked(True)
          self.pushButtonS.setChecked(True)
          self.do_S = 1
#          self.pushButtonS0.setChecked(True)
          self.pushButtonS0.setChecked(True)
          self.do_S0 = 1

        else:

#          self.pushButtonS.setChecked(False)
          self.pushButtonS.setChecked(False)
          self.do_S = 0
          self.pushButtonS0.setChecked(False)
          self.do_S0 = 0

        self.Run()

    def P_clicked(self):
        if ( self.do_P == 0):
          self.DisableAllP()

          self.pushButtonP.setChecked(True)
          self.do_P = 1

          self.pushButtonP_1.setChecked(True)
          self.do_P_1 = 1
          self.pushButtonP0.setChecked(True)
          self.do_P0 = 1
          self.pushButtonP1.setChecked(True)
          self.do_P1 = 1

          # self.select_mode=0
          # self.selected_l = 1
          #self.Run()


        else:

          self.pushButtonP.setChecked(False)
          self.do_P = 0

          self.pushButtonP_1.setChecked(False)
          self.do_P_1 = 0
          self.pushButtonP0.setChecked(False)
          self.do_P0 = 0
          self.pushButtonP1.setChecked(False)
          self.do_P1 = 0

        self.Run()

    def D_clicked(self):
        if ( self.do_D == 0):
          self.DisableAllD()

          self.pushButtonD.setChecked(True)
          self.do_D = 1

          self.pushButtonD_2.setChecked(True)
          self.do_D_2 = 1
          self.pushButtonD_1.setChecked(True)
          self.do_D_1 = 1
          self.pushButtonD0.setChecked(True)
          self.do_D0 = 1
          self.pushButtonD1.setChecked(True)
          self.do_D1 = 1
          self.pushButtonD2.setChecked(True)
          self.do_D2 = 1

          # self.select_mode=0
          # self.selected_l = 2
          #self.Run()


        else:

          self.pushButtonD.setChecked(False)
          self.do_D = 0

          self.pushButtonD_2.setChecked(False)
          self.do_D_2 = 0
          self.pushButtonD_1.setChecked(False)
          self.do_D_1 = 0
          self.pushButtonD0.setChecked(False)
          self.do_D0 = 0
          self.pushButtonD1.setChecked(False)
          self.do_D1 = 0
          self.pushButtonD2.setChecked(False)
          self.do_D2 = 0

        self.Run()

    def F_clicked(self):
        if ( self.do_F == 0):
          self.DisableAllF()

          self.pushButtonF.setChecked(True)
          self.do_F = 1

          self.pushButtonF_3.setChecked(True)
          self.do_F_3 = 1
          self.pushButtonF_2.setChecked(True)
          self.do_F_2 = 1
          self.pushButtonF_1.setChecked(True)
          self.do_F_1 = 1
          self.pushButtonF0.setChecked(True)
          self.do_F0 = 1
          self.pushButtonF1.setChecked(True)
          self.do_F1 = 1
          self.pushButtonF2.setChecked(True)
          self.do_F2 = 1
          self.pushButtonF3.setChecked(True)
          self.do_F3 = 1

          # self.select_mode=0
          # self.selected_l = 3


        else:

          self.pushButtonF.setChecked(False)
          self.do_F = 0

          self.pushButtonF_3.setChecked(False)
          self.do_F_3 = 0
          self.pushButtonF_2.setChecked(False)
          self.do_F_2 = 0
          self.pushButtonF_1.setChecked(False)
          self.do_F_1 = 0
          self.pushButtonF0.setChecked(False)
          self.do_F0 = 0
          self.pushButtonF1.setChecked(False)
          self.do_F1 = 0
          self.pushButtonF2.setChecked(False)
          self.do_F2 = 0
          self.pushButtonF3.setChecked(False)
          self.do_F3 = 0

# m buttons

    def S0_clicked(self):
        if ( self.do_S0 == 0):

          self.pushButtonS0.setChecked(True)
          self.do_S0 = 1

          self.pushButtonS.setChecked(True)

        else:

          self.pushButtonS0.setChecked(False)
          self.do_S0 = 0

        self.Run()
    def P_1_clicked(self):
        if ( self.do_P_1 == 0):

          self.pushButtonP_1.setChecked(True)
          self.do_P_1 = 1

          self.pushButtonP.setChecked(True)

        else:

          self.pushButtonP_1.setChecked(False)
          self.do_P_1 = 0

          if(self.do_P_1==0 and self.do_P0==0 and self.do_P1==0):
            self.pushButtonP.setChecked(False)

        self.Run()

    def P0_clicked(self):
        if ( self.do_P0 == 0):

          self.pushButtonP0.setChecked(True)
          self.do_P0 = 1

          self.pushButtonP.setChecked(True)

        else:

          self.pushButtonP0.setChecked(False)
          self.do_P0 = 0

          if(self.do_P_1==0 and self.do_P0==0 and self.do_P1==0):
            self.pushButtonP.setChecked(False)

        self.Run()


    def P1_clicked(self):
        if ( self.do_P1 == 0):

          self.pushButtonP1.setChecked(True)
          self.do_P1 = 1

          self.pushButtonP.setChecked(True)

        else:

          self.pushButtonP1.setChecked(False)
          self.do_P1 = 0

          if(self.do_P_1==0 and self.do_P0==0 and self.do_P1==0):
            self.pushButtonP.setChecked(False)

        self.Run()

    def D_2_clicked(self):
        if ( self.do_D_2 == 0):

          self.pushButtonD_2.setChecked(True)
          self.do_D_2 = 1

          self.pushButtonD.setChecked(True)

        else:

          self.pushButtonD_2.setChecked(False)
          self.do_D_2 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setChecked(False)

        self.Run()


    def D_1_clicked(self):
        if ( self.do_D_1 == 0):

          self.pushButtonD_1.setChecked(True)
          self.do_D_1 = 1

          self.pushButtonD.setChecked(True)

        else:

          self.pushButtonD_1.setChecked(False)
          self.do_D_1 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setChecked(False)

        self.Run()

    def D0_clicked(self):
        if ( self.do_D0 == 0):

          self.pushButtonD0.setChecked(True)
          self.do_D0 = 1

          self.pushButtonD.setChecked(True)

        else:

          self.pushButtonD0.setChecked(False)
          self.do_D0 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setChecked(False)

        self.Run()

    def D1_clicked(self):
        if ( self.do_D1 == 0):

          self.pushButtonD1.setChecked(True)
          self.do_D1 = 1

          self.pushButtonD.setChecked(True)

        else:

          self.pushButtonD1.setChecked(False)
          self.do_D1 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setChecked(False)

        self.Run()

    def D2_clicked(self):
        if ( self.do_D2 == 0):

          self.pushButtonD2.setChecked(True)
          self.do_D2 = 1

          self.pushButtonD.setChecked(True)


        else:

          self.pushButtonD2.setChecked(False)
          self.do_D2 = 0

          if(self.do_D_2==0 and self.do_D_1==0 and self.do_D0==0 and self.do_D1==0 and self.do_D2==0):
            self.pushButtonD.setChecked(False)

        self.Run()

    def F_3_clicked(self):
        if ( self.do_F_3 == 0):

          self.pushButtonF_3.setChecked(True)
          self.do_F_3 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF_3.setChecked(False)
          self.do_F_3 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)


    def F_2_clicked(self):
        if ( self.do_F_2 == 0):

          self.pushButtonF_2.setChecked(True)
          self.do_F_2 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF_2.setChecked(False)
          self.do_F_2 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)


    def F_1_clicked(self):
        if ( self.do_F_1 == 0):

          self.pushButtonF_1.setChecked(True)
          self.do_F_1 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF_1.setChecked(False)
          self.do_F_1 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)


    def F0_clicked(self):
        if ( self.do_F0 == 0):

          self.pushButtonF0.setChecked(True)
          self.do_F0 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF0.setChecked(False)
          self.do_F0 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)


    def F1_clicked(self):
        if ( self.do_F1 == 0):

          self.pushButtonF1.setChecked(True)
          self.do_F1 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF1.setChecked(False)
          self.do_F1 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)


    def F2_clicked(self):
        if ( self.do_F2 == 0):

          self.pushButtonF2.setChecked(True)
          self.do_F2 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF2.setChecked(False)
          self.do_F2 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)


    def F3_clicked(self):
        if ( self.do_F3 == 0):

          self.pushButtonF3.setChecked(True)
          self.do_F3 = 1

          self.pushButtonF.setChecked(True)

        else:

          self.pushButtonF3.setChecked(False)
          self.do_F3 = 0

          if(self.do_F_3==0 and self.do_F_2==0 and self.do_F_1==0 and self.do_F0==0 and self.do_F1==0 and self.do_F2==0 and self.do_F3==0):
            self.pushButtonF.setChecked(False)



    def printButtonPressed(self):
        # This is executed when the button is pressed
        print('printButtonPressed')

app = QtWidgets.QApplication(sys.argv)
window = Ui()

#plt.show()

app.exec_()
